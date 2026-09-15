"""Cluster genomes into pangenomes using skani or nucmer"""
from __future__ import print_function, division
import sys
import os
import re
import gzip
import shutil
import tarfile
import subprocess
import argparse
from itertools import combinations
from multiprocessing import cpu_count
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from tqdm import tqdm
from loguru import logger
from pyexeggutor import RunShellCommand, format_header
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header


def _check_status_quiet(step, expected_outputs=None):
    """Validate command success without printing to stderr.

    Replaces step.check_status() which unconditionally prints
    'Command Successful: ...' to stderr for every command.
    """
    if step.returncode_ != 0:
        raise subprocess.CalledProcessError(
            returncode=step.returncode_,
            cmd=step.command,
        )
    if expected_outputs:
        for filepath in expected_outputs:
            if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
                raise FileNotFoundError(f"Expected output not found or empty: {filepath}")


def archive_nucmer_results(nucmer_work_dir, archive_dir):
    """Create .tar.gz archives of nucmer intermediate files grouped by type."""
    os.makedirs(archive_dir, exist_ok=True)

    file_groups = {
        "delta": [".delta"],
        "filtered_delta": [".1delta", ".mdelta"],
        "reports": [".report"],
        "coords": [".1coords", ".mcoords"],
        "snps": [".snps"],
        "diff": [".rdiff", ".qdiff", ".unref", ".unqry"],
    }

    for group_name, extensions in file_groups.items():
        matching_files = []
        for ext in extensions:
            matching_files.extend(
                f for f in os.listdir(nucmer_work_dir)
                if f.endswith(ext)
            )
        if not matching_files:
            continue

        archive_path = os.path.join(archive_dir, f"nucmer_results.{group_name}.tar.gz")
        with tarfile.open(archive_path, "w:gz") as tar:
            for fname in sorted(matching_files):
                tar.add(
                    os.path.join(nucmer_work_dir, fname),
                    arcname=fname,
                )
        logger.info(f"Archived {len(matching_files)} files → {os.path.basename(archive_path)}")


def get_basename_from_filepath(filepath):
    """Replicate the basename extraction logic from edgelist-to-clusters.py --basename.

    This must stay in sync with the get_basename() function in edgelist-to-clusters.py
    (lines 102-106) to correctly predict what IDs will appear in the processed edge list.
    """
    _, fn = os.path.split(filepath)
    if fn.endswith(".gz"):
        fn = fn[:-3]
    return ".".join(fn.split(".")[:-1])


def get_file_extension(filepath):
    """Extract extension from filepath, preserving compound extensions like .fa.gz"""
    basename = os.path.basename(filepath)
    if basename.endswith(".gz"):
        inner = basename[:-3]
        parts = inner.split(".")
        if len(parts) > 1:
            return "." + parts[-1] + ".gz"
        return ".gz"
    else:
        parts = basename.split(".")
        if len(parts) > 1:
            return "." + parts[-1]
        return ""


def create_genome_symlinks(genome_id_to_filepath, tmp_directory):
    """Create symlinks named by genome_id when IDs don't match filename basenames.

    This ensures that skani output (which uses filepaths) will produce basenames
    matching the manifest's genome IDs when processed by edgelist-to-clusters.py --basename.

    Args:
        genome_id_to_filepath: OrderedDict mapping genome_id -> filepath
        tmp_directory: Path to tmp directory for symlink storage

    Returns:
        tuple: (genome_id_to_filepath_for_skani, used_symlinks)
            - genome_id_to_filepath_for_skani: OrderedDict with symlink paths (or original if no symlinks needed)
            - used_symlinks: bool indicating whether symlinks were created
    """
    # Check if any genome ID differs from what --basename would derive
    needs_symlinks = False
    for genome_id, filepath in genome_id_to_filepath.items():
        derived_id = get_basename_from_filepath(filepath)
        if genome_id != derived_id:
            needs_symlinks = True
            break

    if not needs_symlinks:
        return genome_id_to_filepath, False

    # Create symlink directory
    symlink_dir = os.path.join(tmp_directory, "genome_symlinks")
    os.makedirs(symlink_dir, exist_ok=True)

    # Create symlinks
    genome_id_to_symlink = OrderedDict()
    for genome_id, filepath in genome_id_to_filepath.items():
        abs_target = os.path.abspath(filepath)
        ext = get_file_extension(filepath)
        symlink_name = genome_id + ext
        symlink_path = os.path.join(symlink_dir, symlink_name)

        # Remove existing symlink if present (idempotent)
        if os.path.lexists(symlink_path):
            os.remove(symlink_path)

        os.symlink(abs_target, symlink_path)
        genome_id_to_symlink[genome_id] = symlink_path

    return genome_id_to_symlink, True


def parse_input(input_path, genome_extension=None):
    """Parse input as either manifest or simple list, autodetecting format
    
    Accepted organism_types: {prokaryotic, eukaryotic, viral}

    Supports:
    - Manifest (3 cols): [organism_type, id_genome, genome_filepath]
    - Manifest (5+ cols): [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, ...]
    - Simple list: one genome filepath per line
    """
    genome_id_to_filepath = OrderedDict()
    
    # Handle stdin
    if input_path == "stdin" or input_path == "-":
        input_handle = sys.stdin
    else:
        if not os.path.exists(input_path):
            raise FileNotFoundError(f"Input file not found: {input_path}")
        input_handle = open(input_path, 'r')
    
    try:
        # Read first line to detect format
        first_line = input_handle.readline().strip()
        
        # Check if it's tab-separated (manifest) or not (simple list)
        if '\t' in first_line:
            # It's a manifest - parse with pandas
            if input_path == "stdin" or input_path == "-":
                # Can't seek stdin, so we need to reconstruct
                from io import StringIO
                content = first_line + '\n' + input_handle.read()
                df = pd.read_csv(StringIO(content), sep="\t", header=None)
            else:
                input_handle.close()
                df = pd.read_csv(input_path, sep="\t", header=None)
            
            if df.shape[1] == 3:
                # Batch mode (simple): [organism_type, id_genome, genome_filepath]
                df.columns = ["organism_type", "id_genome", "genome_filepath"]
                for _, row in df.iterrows():
                    genome_id_to_filepath[row["id_genome"]] = row["genome_filepath"]
            
            elif df.shape[1] == 4:
                # Batch mode (full): [organism_type, id_genome, genome_filepath, protein_filepath]
                df.columns = ["organism_type", "id_genome", "genome_filepath", "protein_filepath"]
                for _, row in df.iterrows():
                    genome_id_to_filepath[row["id_genome"]] = row["genome_filepath"]
                    
            elif df.shape[1] == 5:
                # Batch mode with CDS: [organism_type, id_genome, genome_filepath, protein_filepath, cds_filepath]
                for _, row in df.iterrows():
                    genome_id_to_filepath[row[1]] = row[2]

            elif df.shape[1] >= 6:
                # VEBA mode: [organism_type, id_sample, id_genome, genome_filepath, ...]
                for _, row in df.iterrows():
                    id_genome = row[2]
                    genome_filepath = row[3]
                    genome_id_to_filepath[id_genome] = genome_filepath
            else:
                raise ValueError(f"Manifest must have 3, 4, or ≥5 columns. Got {df.shape[1]}.")
        
        else:
            # It's a simple list - process line by line
            # First pass: collect all filepaths and check extensions
            all_filepaths = [first_line] + [line.strip() for line in input_handle]
            all_filepaths = [fp.strip() for fp in all_filepaths if fp.strip()]
            
            # Validate extension if provided
            if genome_extension:
                ext = genome_extension
                if not ext.startswith("."):
                    ext = "." + ext
                
                # Check how many files have the provided extension
                files_with_ext = [fp for fp in all_filepaths if os.path.basename(fp).endswith(ext)]
                
                if len(files_with_ext) == 0:
                    # No files have this extension - this is an error
                    actual_extensions = set(os.path.splitext(os.path.basename(fp))[1] for fp in all_filepaths)
                    raise ValueError(
                        f"Extension mismatch: You provided extension '{genome_extension}' but none of the "
                        f"{len(all_filepaths)} input files have this extension.\n"
                        f"Actual extensions found: {', '.join(sorted(actual_extensions))}\n"
                        f"Hint: Use '-x {list(actual_extensions)[0][1:]}' for '{list(actual_extensions)[0]}' files"
                    )
                elif len(files_with_ext) < len(all_filepaths):
                    # Some files don't have this extension - warning
                    from loguru import logger
                    logger.warning(
                        f"Extension mismatch: {len(files_with_ext)}/{len(all_filepaths)} files have "
                        f"extension '{ext}'. Files without this extension will be parsed differently."
                    )
            
            # Second pass: process filepaths
            for filepath in all_filepaths:
                basename = os.path.basename(filepath)
                
                if genome_extension:
                    ext = genome_extension
                    if not ext.startswith("."):
                        ext = "." + ext
                    genome_id = basename[:-len(ext)] if basename.endswith(ext) else basename
                else:
                    genome_id = get_basename_from_filepath(filepath)
                
                genome_id_to_filepath[genome_id] = filepath
    
    finally:
        if input_handle != sys.stdin:
            input_handle.close()
    
    return genome_id_to_filepath

def decompress_genomes(genome_id_to_filepath, tmp_directory):
    """Decompress .gz genome files to a temporary directory.

    Returns an updated genome_id -> filepath mapping where .gz files
    are replaced with paths to their decompressed copies. Non-gzipped
    files are returned as-is.
    """
    genome_id_to_decompressed = OrderedDict()
    needs_decompress = any(fp.endswith(".gz") for fp in genome_id_to_filepath.values())

    if not needs_decompress:
        return genome_id_to_filepath

    decompressed_dir = os.path.join(tmp_directory, "decompressed_genomes")
    os.makedirs(decompressed_dir, exist_ok=True)
    logger.info(f"Decompressing .gz genomes to {decompressed_dir}")

    for genome_id, filepath in genome_id_to_filepath.items():
        if filepath.endswith(".gz"):
            ext = get_file_extension(filepath)
            if ext.endswith(".gz"):
                ext = ext[:-3]
            decompressed_path = os.path.join(decompressed_dir, genome_id + ext)
            with gzip.open(filepath, 'rb') as f_in, open(decompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            genome_id_to_decompressed[genome_id] = decompressed_path
        else:
            genome_id_to_decompressed[genome_id] = filepath

    return genome_id_to_decompressed


def parse_dnadiff_report(report_path, ref_id, qry_id, identity_type="1-to-1"):
    """Parse a dnadiff .report file to extract ANI and alignment fractions.

    Returns (ref_id, qry_id, ani, af_ref, af_query) or None if parsing fails.
    """
    af_ref = None
    af_query = None
    ani = None
    current_section = None

    with open(report_path, 'r') as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("1-to-1"):
                current_section = "1-to-1"
                continue
            elif stripped.startswith("M-to-M"):
                current_section = "M-to-M"
                continue

            if stripped.startswith("AlignedBases"):
                pcts = re.findall(r'\((\d+\.?\d*)%\)', stripped)
                if len(pcts) >= 2:
                    af_ref = float(pcts[0])
                    af_query = float(pcts[1])

            if stripped.startswith("AvgIdentity") and current_section == identity_type:
                parts = stripped.split()
                if len(parts) >= 2:
                    ani = float(parts[1])

    if ani is not None and af_ref is not None and af_query is not None:
        return (ref_id, qry_id, ani, af_ref, af_query)

    logger.warning(f"Could not fully parse dnadiff report: {report_path}")
    return None


def run_nucmer_pair(ref_id, qry_id, ref_path, qry_path, work_dir,
                    nucmer_options="", identity_type="1-to-1", n_threads=1):
    """Run nucmer + dnadiff for a single genome pair.

    Returns (ref_id, qry_id, ani, af_ref, af_query) or None.
    """
    pair_prefix = f"{ref_id}__vs__{qry_id}"
    prefix = os.path.join(work_dir, pair_prefix)

    cmd = ["nucmer", "-p", prefix, "-t", str(n_threads)]
    if nucmer_options:
        cmd.extend(nucmer_options.split())
    cmd.extend([ref_path, qry_path])

    delta_file = f"{prefix}.delta"
    report_file = f"{prefix}.report"

    step = RunShellCommand(
        command=cmd,
        name=f"nucmer:{pair_prefix}",
    ).run()
    _check_status_quiet(step, expected_outputs=[delta_file])

    cmd_dnadiff = ["dnadiff", "-d", delta_file, "-p", prefix]
    step = RunShellCommand(
        command=cmd_dnadiff,
        name=f"dnadiff:{pair_prefix}",
    ).run()
    _check_status_quiet(step, expected_outputs=[report_file])

    return parse_dnadiff_report(report_file, ref_id, qry_id, identity_type)


def run_nucmer_all_vs_all(genome_id_to_filepath, work_dir, output_edgelist,
                          nucmer_options="", identity_type="1-to-1",
                          n_threads_per_task=1, n_concurrent_tasks=1):
    """Run nucmer + dnadiff for all pairwise genome comparisons.

    Writes a 5-column edge list: [id_ref, id_qry, ANI, AF_ref, AF_query].
    Returns the DataFrame of results.
    """
    pairs = list(combinations(genome_id_to_filepath.keys(), 2))
    n_pairs = len(pairs)
    logger.info(f"Running {n_pairs} pairwise nucmer comparisons ({len(genome_id_to_filepath)} genomes)")

    results = []

    def _run_pair(ref_id, qry_id):
        return run_nucmer_pair(
            ref_id, qry_id,
            genome_id_to_filepath[ref_id],
            genome_id_to_filepath[qry_id],
            work_dir, nucmer_options, identity_type, n_threads_per_task
        )

    if n_concurrent_tasks > 1:
        with ThreadPoolExecutor(max_workers=n_concurrent_tasks) as executor:
            futures = {
                executor.submit(_run_pair, ref_id, qry_id): (ref_id, qry_id)
                for ref_id, qry_id in pairs
            }
            for future in tqdm(as_completed(futures), total=n_pairs,
                               desc="Nucmer pairwise", unit=" pairs"):
                result = future.result()
                if result is not None:
                    results.append(result)
    else:
        for ref_id, qry_id in tqdm(pairs, desc="Nucmer pairwise", unit=" pairs"):
            result = _run_pair(ref_id, qry_id)
            if result is not None:
                results.append(result)

    df_edges = pd.DataFrame(results, columns=["id_ref", "id_qry", "ANI", "AF_ref", "AF_query"])
    df_edges.to_csv(output_edgelist, sep="\t", index=False, header=False)
    logger.info(f"Wrote {len(results)} edges to {output_edgelist}")

    return df_edges


def archive_dotplot_auxiliary_files(dotplot_dir):
    """Archive mummerplot auxiliary files (.gp, .fplot, .rplot) and remove originals."""
    aux_extensions = [".gp", ".fplot", ".rplot"]
    matching_files = []
    for ext in aux_extensions:
        matching_files.extend(
            f for f in os.listdir(dotplot_dir)
            if f.endswith(ext)
        )

    if not matching_files:
        return

    archive_path = os.path.join(dotplot_dir, "dotplot_auxiliary_files.tar.gz")
    with tarfile.open(archive_path, "w:gz") as tar:
        for fname in sorted(matching_files):
            tar.add(
                os.path.join(dotplot_dir, fname),
                arcname=fname,
            )

    for fname in matching_files:
        os.remove(os.path.join(dotplot_dir, fname))

    logger.info(f"Archived {len(matching_files)} auxiliary files → {os.path.basename(archive_path)}")


def generate_dotplots(passing_pairs, genome_id_to_decompressed, nucmer_work_dir,
                      dotplot_dir, dotplot_format="pdf", nucmer_options="", n_threads=1):
    """Generate dot plots for threshold-passing genome pairs.

    For nucmer backend, delta files already exist in nucmer_work_dir.
    For skani backend, runs nucmer first to produce delta files.
    """
    os.makedirs(dotplot_dir, exist_ok=True)

    for ref_id, qry_id in tqdm(passing_pairs, desc="Generating dot plots", unit=" plots"):
        pair_prefix = f"{ref_id}__vs__{qry_id}"
        delta_file = os.path.join(nucmer_work_dir, f"{pair_prefix}.delta")

        if not os.path.exists(delta_file):
            os.makedirs(nucmer_work_dir, exist_ok=True)
            cmd = ["nucmer", "-p", os.path.join(nucmer_work_dir, pair_prefix),
                   "-t", str(n_threads)]
            if nucmer_options:
                cmd.extend(nucmer_options.split())
            cmd.extend([
                genome_id_to_decompressed[ref_id],
                genome_id_to_decompressed[qry_id]
            ])
            step = RunShellCommand(
                command=cmd,
                name=f"nucmer_dotplot:{pair_prefix}",
            ).run()
            _check_status_quiet(step, expected_outputs=[delta_file])

        dotplot_prefix = os.path.join(dotplot_dir, pair_prefix)
        cmd = [
            "mummerplot", delta_file,
            "-p", dotplot_prefix,
            "-t", dotplot_format,
            "--large",
        ]
        step = RunShellCommand(
            command=cmd,
            name=f"mummerplot:{pair_prefix}",
        ).run()
        _check_status_quiet(step)

    logger.info(f"Generated {len(passing_pairs)} dot plots in {dotplot_dir}")


def register_parser(subparsers):
    """Register the cluster-genomes subcommand parser"""
    parser = subparsers.add_parser(
        'cluster-genomes',
        help='Cluster genomes into pangenomes using skani or nucmer',
        description='Cluster genomes into pangenomes based on Average Nucleotide Identity (ANI)\n'
                    'Supports two backends: skani (fast, default) and nucmer (MUMmer4, pairwise)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument(
        "-i", "--input",
        type=str,
        default="stdin",
        help="Input: manifest file OR list of genome filepaths (one per line) [Default: stdin]"
    )
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenomium_output/genome_clustering", help="Output directory [Default: pangenomium_output/genome_clustering]")
    parser_io.add_argument("-x", "--genome_extension", type=str, help="Genome extension for parsing IDs from list (e.g., 'fa.gz')")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--n_threads", type=int, default=1, help="Number of threads [Default: 1]")
    parser_utility.add_argument("--keep_temporary", action="store_true", help="Keep temporary directories (default: remove after completion)")

    # Algorithm selection
    parser_algorithm = parser.add_argument_group('Algorithm selection')
    parser_algorithm.add_argument("--genome_clustering_algorithm", type=str, default="skani",
        choices=["skani", "nucmer"],
        help="Algorithm for genome clustering:\n"
             "  skani:   Fast ANI via skani triangle (default)\n"
             "  nucmer:  ANI via MUMmer4 nucmer + dnadiff (pairwise)\n"
             "[Default: skani]")

    # ANI threshold arguments (shared by both backends)
    parser_ani = parser.add_argument_group('ANI threshold arguments')
    parser_ani.add_argument("--ani_threshold", type=float, default=95.0, help="ANI threshold [Default: 95.0]")
    parser_ani.add_argument("--minimum_af", type=float, default=50.0, help="Minimum alignment fraction [Default: 50.0]")
    parser_ani.add_argument("--af_mode", type=str, default="relaxed", choices=["relaxed", "strict"], help="AF mode [Default: relaxed]")

    # Skani arguments
    parser_skani = parser.add_argument_group('Skani arguments (only with --genome_clustering_algorithm skani)')
    parser_skani.add_argument("--skani_preset", type=str, help="Skani preset")
    parser_skani.add_argument("--skani_options", type=str, default="", help="Additional skani options")

    # Nucmer arguments
    parser_nucmer = parser.add_argument_group('Nucmer arguments (only with --genome_clustering_algorithm nucmer)')
    parser_nucmer.add_argument("--nucmer_options", type=str, default="",
        help="Additional nucmer options (e.g., '--maxmatch')")
    parser_nucmer.add_argument("--n_concurrent_nucmer_tasks", type=int, default=1,
        help="Number of concurrent nucmer pairwise comparisons [Default: 1]")
    parser_nucmer.add_argument("--nucmer_identity_type", type=str, default="1-to-1",
        choices=["1-to-1", "M-to-M"],
        help="Which dnadiff AvgIdentity to use as ANI:\n"
             "  1-to-1: Best bidirectional alignment identity (conservative)\n"
             "  M-to-M: Many-to-many alignment identity (more permissive)\n"
             "[Default: 1-to-1]")

    # Dot plot arguments (both backends)
    parser_dotplot = parser.add_argument_group('Dot plot arguments')
    parser_dotplot.add_argument("--generate_dotplots", action="store_true",
        help="Generate mummerplot dot plots for threshold-passing genome pairs")
    parser_dotplot.add_argument("--dotplot_format", type=str, default="pdf",
        choices=["pdf", "png", "ps", "svg"],
        help="Dot plot output format [Default: pdf]")
    
    # Clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--organism_type", type=str, choices=["prokaryotic", "eukaryotic", "viral"], help="Organism type (required with --prepend_organism_code)")
    parser_clustering.add_argument("--prepend_organism_code", action="store_true", help="Prepend organism code to cluster prefix (P/E/V for prokaryotic/eukaryotic/viral)")
    parser_clustering.add_argument("--cluster_prefix", type=str, default="SLC-", help="Cluster prefix [Default: 'SLC-']")
    parser_clustering.add_argument("--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")
    parser_clustering.add_argument("--cluster_prefix_zfill", type=int, default=0, help="Prefix zfill [Default: 0]")
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    parser_clustering.add_argument("--identifiers", type=str, help="Identifiers file")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def run(args):
    """Execute cluster-genomes command"""

    if args.n_threads == -1:
        args.n_threads = cpu_count()

    # Validate organism code usage
    if args.prepend_organism_code and not args.organism_type:
        logger.error("--prepend_organism_code requires --organism_type to be specified")
        return 1

    # Prepend organism code to cluster prefix if requested
    cluster_prefix = args.cluster_prefix
    if args.prepend_organism_code and args.organism_type:
        organism_code = args.organism_type[0].upper()  # P, E, or V
        cluster_prefix = organism_code + cluster_prefix
        logger.info(f"Prepending organism code '{organism_code}' to cluster prefix: {cluster_prefix}")

    # Setup directories with command-specific subdirectory
    directories = setup_directories(args.output_directory, subdirectory="genome_clustering")

    # Create additional output subdirectories
    os.makedirs(os.path.join(directories["output"], "serialization"), exist_ok=True)
    os.makedirs(os.path.join(directories["output"], "representatives"), exist_ok=True)

    # Setup logger
    setup_logger(directories["log"], "cluster_genomes.log")

    # Print info
    logger.info("="*80)
    logger.info("pangenomium cluster-genomes")
    logger.info("="*80)
    print_header(
        version=__version__,
        n_jobs=args.n_threads,
        additional_info={
            "Algorithm": args.genome_clustering_algorithm,
            "ANI threshold": args.ani_threshold,
            "Minimum AF": args.minimum_af,
        }
    )

    # Parse input
    logger.info("-"*80)
    logger.info("Parsing input")
    logger.info("-"*80)
    input_source = "stdin" if args.input in ["stdin", "-"] else args.input
    logger.info(f"Input source: {input_source}")

    genome_id_to_filepath = parse_input(args.input, args.genome_extension)
    logger.info(f"Genomes: {len(genome_id_to_filepath)}")

    # Write genome identifiers list (shared by both backends, used by edgelist-to-clusters)
    genome_identifiers_filepath = os.path.join(directories["intermediate"], "genome_identifiers.list")
    with open(genome_identifiers_filepath, "w") as f:
        for genome_id in genome_id_to_filepath.keys():
            print(genome_id, file=f)

    logger.info("")

    # Paths used by both branches
    ani_edgelist_processed = os.path.join(directories["intermediate"], "ani_edgelist_processed.tsv")
    nucmer_work_dir = os.path.join(directories["intermediate"], "nucmer_results")

    # ==========================================
    # Step 1: Compute all-vs-all ANI
    # ==========================================
    if args.genome_clustering_algorithm == "skani":
        # --- Skani backend ---
        logger.info("="*80)
        logger.info("Step 1: Running skani triangle")
        logger.info("="*80)

        # Create symlinks if genome IDs don't match filename basenames
        genome_id_to_filepath_for_skani, used_symlinks = create_genome_symlinks(
            genome_id_to_filepath, directories["tmp"]
        )
        if used_symlinks:
            logger.info("Created genome symlinks (genome IDs differ from filenames)")

        # Write genome list (filepaths for skani)
        genome_list_filepath = os.path.join(directories["intermediate"], "genome_list.txt")
        with open(genome_list_filepath, "w") as f:
            for genome_id, filepath in genome_id_to_filepath_for_skani.items():
                print(filepath, file=f)

        ani_edgelist = os.path.join(directories["intermediate"], "skani-triangle_results.tsv")

        cmd = [
            "skani", "triangle",
            "-l", genome_list_filepath,
            "-E",
            "-t", str(args.n_threads),
            "-o", ani_edgelist,
        ]

        if args.skani_preset:
            cmd.extend(["--preset", args.skani_preset])

        if args.skani_options:
            cmd.append(args.skani_options)

        step = RunShellCommand(
            command=cmd,
            name="skani",
            validate_output_filepaths=[ani_edgelist]
        ).run()
        step.check_status()
        logger.info("")

        # Preprocess skani output: extract first 5 columns and remove header
        logger.info("Preprocessing skani output (extracting columns 1-5)")
        df_ani = pd.read_csv(ani_edgelist, sep="\t", usecols=[0,1,2,3,4])
        df_ani.to_csv(ani_edgelist_processed, sep="\t", index=False, header=False)
        logger.info(f"Processed: {df_ani.shape[0]} edges")
        logger.info("")

    elif args.genome_clustering_algorithm == "nucmer":
        # --- Nucmer backend ---
        logger.info("="*80)
        logger.info("Step 1: Running nucmer + dnadiff pairwise comparisons")
        logger.info("="*80)

        os.makedirs(nucmer_work_dir, exist_ok=True)

        # Decompress .gz genomes (nucmer cannot read gzipped FASTA)
        genome_id_to_decompressed = decompress_genomes(
            genome_id_to_filepath, directories["tmp"]
        )

        run_nucmer_all_vs_all(
            genome_id_to_filepath=genome_id_to_decompressed,
            work_dir=nucmer_work_dir,
            output_edgelist=ani_edgelist_processed,
            nucmer_options=args.nucmer_options,
            identity_type=args.nucmer_identity_type,
            n_threads_per_task=args.n_threads,
            n_concurrent_tasks=args.n_concurrent_nucmer_tasks,
        )

        # Archive nucmer intermediate files by type
        logger.info("Archiving nucmer intermediate files")
        archive_nucmer_results(nucmer_work_dir, os.path.join(directories["intermediate"], "archives"))
        logger.info("")

    # ========================
    # Step 2: Compile clusters
    # ========================
    logger.info("="*80)
    logger.info("Step 2: Compiling genome clusters")
    logger.info("="*80)

    genome_clusters = os.path.join(directories["output"], "genomes_to_pangenomes.tsv.gz")

    cmd = [
        "edgelist-to-clusters.py",
        "-i", ani_edgelist_processed,
        "-o", genome_clusters,
        "-t", str(args.ani_threshold),
        "-a", str(args.minimum_af),
        "-m", args.af_mode,
        "--cluster_prefix", cluster_prefix,
        "--cluster_prefix_zfill", str(args.cluster_prefix_zfill),
        "--cluster_label_mode", args.cluster_label_mode,
        "--identifiers", genome_identifiers_filepath,
        "-g", os.path.join(directories["output"], "serialization", "genome_clusters.graph.pkl.gz"),
        "-d", os.path.join(directories["output"], "serialization", "genome_clusters.dict.pkl.gz"),
        "-r", os.path.join(directories["output"], "representatives", "genome_representatives.tsv.gz"),
    ]

    # skani outputs filepaths in the edge list; nucmer outputs genome IDs directly
    if args.genome_clustering_algorithm == "skani":
        cmd.append("--basename")

    if args.cluster_suffix:
        cmd.extend(["--cluster_suffix", args.cluster_suffix])

    if args.no_singletons:
        cmd.append("--no_singletons")

    if args.identifiers:
        logger.warning("Ignoring --identifiers argument, using auto-generated identifiers list")

    step = RunShellCommand(
        command=cmd,
        name="compile",
        validate_output_filepaths=[genome_clusters]
    ).run()
    step.check_status()
    logger.info("")

    # ==============================
    # Step 3: Dot plots (optional)
    # ==============================
    if args.generate_dotplots:
        logger.info("="*80)
        logger.info("Step 3: Generating dot plots")
        logger.info("="*80)

        # Read the edge list and filter by thresholds
        df_edges = pd.read_csv(ani_edgelist_processed, sep="\t", header=None,
                               names=["id_1", "id_2", "ANI", "AF_ref", "AF_query"])

        # For skani, the first two columns are filepaths — convert to genome IDs
        if args.genome_clustering_algorithm == "skani":
            df_edges["id_1"] = df_edges["id_1"].apply(get_basename_from_filepath)
            df_edges["id_2"] = df_edges["id_2"].apply(get_basename_from_filepath)

        # Apply threshold filtering (same logic as edgelist-to-clusters.py)
        mask_ani = df_edges["ANI"] >= args.ani_threshold
        if args.af_mode == "relaxed":
            mask_af = df_edges[["AF_ref", "AF_query"]].max(axis=1) >= args.minimum_af
        else:
            mask_af = (df_edges["AF_ref"] >= args.minimum_af) & (df_edges["AF_query"] >= args.minimum_af)

        df_passing = df_edges[mask_ani & mask_af]
        passing_pairs = list(zip(df_passing["id_1"], df_passing["id_2"]))

        if len(passing_pairs) > 0:
            # Decompress genomes if needed (for nucmer to generate delta files)
            genome_id_to_decompressed = decompress_genomes(
                genome_id_to_filepath, directories["tmp"]
            )

            dotplot_dir = os.path.join(directories["output"], "dotplots")
            generate_dotplots(
                passing_pairs=passing_pairs,
                genome_id_to_decompressed=genome_id_to_decompressed,
                nucmer_work_dir=nucmer_work_dir,
                dotplot_dir=dotplot_dir,
                dotplot_format=args.dotplot_format,
                nucmer_options=getattr(args, 'nucmer_options', ''),
                n_threads=args.n_threads,
            )

            # Archive auxiliary files (.gp, .fplot, .rplot) from dotplot directory
            archive_dotplot_auxiliary_files(dotplot_dir)
        else:
            logger.info("No pairs passed thresholds, skipping dot plot generation")
        logger.info("")

    # ==============================
    # Cleanup temporary directories
    # ==============================
    if not args.keep_temporary and os.path.exists(directories["tmp"]):
        logger.info("Removing temporary directory: {}".format(directories["tmp"]))
        shutil.rmtree(directories["tmp"], ignore_errors=True)

    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    return 0
