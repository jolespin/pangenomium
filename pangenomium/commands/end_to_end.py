"""End-to-end workflow for clustering genomes and proteins with comprehensive output"""
from __future__ import print_function, division
import sys
import os
import argparse
import gzip
from collections import defaultdict, OrderedDict
from multiprocessing import cpu_count
import numpy as np
import pandas as pd
import pyfastx
from loguru import logger
from pyexeggutor import RunShellCommand, format_header, open_file_writer
from tqdm import tqdm
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header

def parse_manifest(input_path, mode="batch"):
    """Parse input manifest
    
    Accepted organism_types: {prokaryotic, eukaryotic, viral}
    
    Batch mode (3-5 columns):
        [organism_type, id_genome, genome_filepath]
        [organism_type, id_genome, genome_filepath, protein_filepath]
        [organism_type, id_genome, genome_filepath, protein_filepath, cds_filepath]
    
    VEBA mode (6+ columns):
        [organism_type, id_sample, id_genome, genome, proteins, cds, ...]
    """
    genome_data = OrderedDict()
    
    # Handle stdin
    if input_path in ["stdin", "-"]:
        input_handle = sys.stdin
    else:
        input_handle = open(input_path, 'r')
    
    try:
        df = pd.read_csv(input_handle, sep="\t", header=None)
        
        if mode == "veba":
            # VEBA mode: 6+ columns
            assert df.shape[1] >= 6, f"VEBA mode requires ≥6 columns, got {df.shape[1]}"
            df = df.iloc[:, :6]
            df.columns = ["organism_type", "id_sample", "id_genome", "genome", "proteins", "cds"]
            
            for _, row in df.iterrows():
                genome_data[row["id_genome"]] = {
                    "organism_type": row["organism_type"],
                    "id_sample": row["id_sample"],
                    "genome": row["genome"],
                    "proteins": row["proteins"],
                    "cds": row["cds"] if pd.notnull(row["cds"]) else None,
                }
        
        else:  # batch mode
            assert df.shape[1] >= 3, f"Batch mode requires ≥3 columns, got {df.shape[1]}"
            
            if df.shape[1] == 3:
                df.columns = ["organism_type", "id_genome", "genome"]
                for _, row in df.iterrows():
                    genome_data[row["id_genome"]] = {
                        "organism_type": row["organism_type"],
                        "id_sample": row["id_genome"],  # Use genome ID as sample
                        "genome": row["genome"],
                        "proteins": None,
                        "cds": None,
                    }
            
            elif df.shape[1] == 4:
                df.columns = ["organism_type", "id_genome", "genome", "proteins"]
                for _, row in df.iterrows():
                    genome_data[row["id_genome"]] = {
                        "organism_type": row["organism_type"],
                        "id_sample": row["id_genome"],
                        "genome": row["genome"],
                        "proteins": row["proteins"],
                        "cds": None,
                    }
            
            elif df.shape[1] >= 5:
                df = df.iloc[:, :5]
                df.columns = ["organism_type", "id_genome", "genome", "proteins", "cds"]
                for _, row in df.iterrows():
                    genome_data[row["id_genome"]] = {
                        "organism_type": row["organism_type"],
                        "id_sample": row["id_genome"],
                        "genome": row["genome"],
                        "proteins": row["proteins"],
                        "cds": row["cds"] if pd.notnull(row["cds"]) else None,
                    }
    
    finally:
        if input_handle != sys.stdin:
            input_handle.close()
    
    return genome_data


def get_protein_cluster_prevalence(df_input: pd.DataFrame):
    """Generate genome × protein_cluster prevalence matrix"""
    # Read Input
    genomes = sorted(df_input.iloc[:, 0].unique())
    clusters = sorted(df_input.iloc[:, 2].unique())

    # Create array
    A = np.zeros((len(genomes), len(clusters)), dtype=int)

    for _, (id_genome, id_protein, id_cluster) in df_input.iterrows():
        i = genomes.index(id_genome)
        j = clusters.index(id_cluster)
        A[i, j] += 1

    # Create output
    df_output = pd.DataFrame(A, index=genomes, columns=clusters)
    df_output.index.name = "id_genome"
    df_output.columns.name = "id_protein_cluster"

    return df_output


def register_parser(subparsers):
    """Register the end-to-end subcommand parser"""
    parser = subparsers.add_parser(
        'end-to-end',
        help='Run complete pangenomium workflow (genomes + proteins) with comprehensive output',
        description='End-to-end workflow: cluster genomes into pangenomes, cluster proteins within each pangenome, and generate comprehensive identifier mappings and statistics',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Input manifest file (see --mode for format)"
    )
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenomium_output", help="Output directory [Default: pangenomium_output]")
    parser_io.add_argument("-m", "--mode", type=str, default="batch", choices=["batch", "veba"], 
                          help="Input mode:\n"
                               "  batch: [organism_type, id_genome, genome, proteins] (4 cols) OR\n"
                               "         [organism_type, id_genome, genome, proteins, cds] (5 cols)\n"
                               "  veba:  [organism_type, id_sample, id_genome, genome, proteins, cds, ...] (6+ cols)\n"
                               "[Default: batch]")
    parser_io.add_argument("--no_core_sequences", action="store_true", help="Don't write core pangenome sequences")
    parser_io.add_argument("--no_representative_sequences", action="store_true", help="Don't write representative sequences")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--n_threads_skani", type=int, default=1, help="Threads for genome clustering (skani) [Default: 1]")
    parser_utility.add_argument("--n_threads_mmseqs_per_task", type=int, default=1, help="Threads per MMseqs2 task [Default: 1]")
    parser_utility.add_argument("--n_concurrent_mmseqs_tasks", type=int, default=1, help="Number of pangenomes to process in parallel [Default: 1]")
    parser_utility.add_argument("--keep_temporary", action="store_true", help="Keep temporary directories (default: remove after completion)")
    
    # Genome clustering arguments
    parser_genome = parser.add_argument_group('Genome clustering arguments')
    parser_genome.add_argument("--ani_threshold", type=float, default=95.0, help="ANI threshold [Default: 95.0]")
    parser_genome.add_argument("--minimum_af", type=float, default=50.0, help="Minimum AF [Default: 50.0]")
    parser_genome.add_argument("--af_mode", type=str, default="relaxed", choices=["relaxed", "strict"], help="AF mode [Default: relaxed]")
    parser_genome.add_argument("--skani_preset", type=str, help="Skani preset")
    parser_genome.add_argument("--skani_options", type=str, default="", help="Additional skani options")
    parser_genome.add_argument("--organism_type", type=str, choices=["prokaryotic", "eukaryotic", "viral"], help="Organism type (required with --prepend_organism_code)")
    parser_genome.add_argument("--prepend_organism_code", action="store_true", help="Prepend organism code to genome cluster prefix (P/E/V)")
    parser_genome.add_argument("--genome_cluster_prefix", type=str, default="SLC-", help="Genome cluster prefix [Default: 'SLC-']")
    
    # Protein clustering arguments
    parser_protein = parser.add_argument_group('Protein clustering arguments')
    parser_protein.add_argument("-a", "--algorithm", type=str, default="mmseqs-cluster", choices=["mmseqs-cluster", "mmseqs-linclust"], help="Algorithm [Default: mmseqs-cluster]")
    parser_protein.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="Identity threshold [Default: 50.0]")
    parser_protein.add_argument("--minimum_coverage_threshold", type=float, default=0.8, help="Coverage threshold [Default: 0.8]")
    parser_protein.add_argument("--separator", type=str, default="_", help="Separator between genome and protein cluster IDs [Default: '_']")
    parser_protein.add_argument("--protein_cluster_prefix", type=str, default="SSPC-", help="Protein cluster prefix [Default: 'SSPC-']")
    parser_protein.add_argument("--mmseqs2_options", type=str, default="", help="MMseqs2 options")
    
    # Pangenome arguments
    parser_pangenome = parser.add_argument_group('Pangenome arguments')
    parser_pangenome.add_argument("--minimum_core_prevalence", type=float, default=1.0, 
                                  help="Minimum ratio of genomes for a protein cluster to be considered core (0.0, 1.0] [Default: 1.0]")
    
    # General clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser


def run(args):
    """Execute end-to-end command with comprehensive output generation"""
    
    # Handle -1 (use all CPUs)
    if args.n_threads_skani == -1:
        args.n_threads_skani = cpu_count()
    if args.n_threads_mmseqs_per_task == -1:
        args.n_threads_mmseqs_per_task = cpu_count()
    if args.n_concurrent_mmseqs_tasks == -1:
        args.n_concurrent_mmseqs_tasks = cpu_count()
    
    assert 0 < args.minimum_core_prevalence <= 1.0, "--minimum_core_prevalence must be (0.0, 1.0]"
    
    # Setup directories (no subdirectory for end-to-end)
    directories = setup_directories(args.output_directory)
    
    # Create all output subdirectories
    os.makedirs(os.path.join(directories["output"], "pangenome_tables"), exist_ok=True)
    os.makedirs(os.path.join(directories["output"], "serialization"), exist_ok=True)
    os.makedirs(os.path.join(directories["output"], "representatives"), exist_ok=True)
    if not args.no_core_sequences:
        os.makedirs(os.path.join(directories["output"], "pangenome_core_sequences"), exist_ok=True)
    
    # Setup logger
    setup_logger(directories["log"], "end_to_end.log")
    
    # Print info
    logger.info("="*80)
    logger.info("pangenomium end-to-end workflow")
    logger.info("="*80)
    print_header(
        version=__version__,
        n_jobs=f"Genome: {args.n_threads_skani} threads | Protein: {args.n_concurrent_mmseqs_tasks} concurrent × {args.n_threads_mmseqs_per_task} threads",
        additional_info={"Mode": args.mode}
    )
    
    # =======================
    # Parse input manifest
    # =======================
    logger.info("="*80)
    logger.info("Step 1: Parsing input manifest")
    logger.info("="*80)
    
    genome_data = parse_manifest(args.input, mode=args.mode)
    logger.info(f"Genomes: {len(genome_data)}")
    
    # Auto-detect organism type from manifest if prepend_organism_code is used without --organism_type
    if args.prepend_organism_code and not args.organism_type:
        organism_types = set(data["organism_type"] for data in genome_data.values())
        if len(organism_types) == 1:
            args.organism_type = list(organism_types)[0]
            logger.info(f"Auto-detected organism type from manifest: {args.organism_type}")
        elif len(organism_types) > 1:
            logger.error(f"--prepend_organism_code requires --organism_type when manifest contains multiple organism types: {organism_types}")
            logger.error("Please specify --organism_type explicitly or ensure all genomes have the same organism type")
            return 1
        else:
            logger.error("No organism types found in manifest")
            return 1
    
    logger.info("")
    
    # Create genome manifest for clustering
    genome_manifest = os.path.join(directories["intermediate"], "genomes_manifest.tsv")
    with open(genome_manifest, "w") as f:
        for id_genome, data in genome_data.items():
            if args.mode == "veba":
                # VEBA format: [organism_type, id_sample, id_genome, genome, proteins, cds]
                cds_val = data["cds"] if data["cds"] is not None else ""
                print(
                    data["organism_type"], 
                    data["id_sample"],
                    id_genome, 
                    data["genome"],
                    data["proteins"],
                    cds_val,
                    sep="\t", 
                    file=f
                )
            else:
                # Batch format: [organism_type, id_genome, genome, proteins, cds]
                if data["cds"] is not None:
                    print(
                        data["organism_type"], 
                        id_genome, 
                        data["genome"],
                        data["proteins"],
                        data["cds"],
                        sep="\t", 
                        file=f
                    )
                else:
                    print(
                        data["organism_type"], 
                        id_genome, 
                        data["genome"],
                        data["proteins"],
                        sep="\t", 
                        file=f
                    )
    
    # ===============================
    # Step 2: Run genome clustering
    # ===============================
    logger.info("="*80)
    logger.info("Step 2: Genome clustering")
    logger.info("="*80)
    
    genome_clustering_dir = os.path.join(directories["project"], "genome_clustering")
    genome_clusters_file = os.path.join(genome_clustering_dir, "output", "genomes_to_pangenomes.tsv.gz")
    
    cmd = [
        "pangenomium", "cluster-genomes",
        "-i", genome_manifest,
        "-o", genome_clustering_dir,
        "--n_threads", str(args.n_threads_skani),
        "--ani_threshold", str(args.ani_threshold),
        "--minimum_af", str(args.minimum_af),
        "--af_mode", args.af_mode,
        "--cluster_prefix", args.genome_cluster_prefix,
        "--cluster_label_mode", args.cluster_label_mode,
    ]
    
    if args.organism_type:
        cmd.extend(["--organism_type", args.organism_type])
    if args.prepend_organism_code:
        cmd.append("--prepend_organism_code")
    if args.skani_preset:
        cmd.extend(["--skani_preset", args.skani_preset])
    if args.skani_options:
        cmd.extend(["--skani_options", args.skani_options])
    if args.no_singletons:
        cmd.append("--no_singletons")
    
    step = RunShellCommand(
        command=cmd,
        name="cluster_genomes",
        validate_output_filepaths=[genome_clusters_file]
    ).run()
    step.check_status()
    logger.info(f"Genome clustering complete")
    logger.info("")
    
    # Create pangenome manifest
    logger.info("="*80)
    logger.info("Step 3: Compiling pangenome manifest")
    logger.info("="*80)
    
    pangenome_manifest = os.path.join(directories["intermediate"], "pangenomes_manifest.tsv")
    
    cmd = [
        "pangenomium", "compile-pangenomes-table",
        "-i", genome_manifest,
        "-c", genome_clusters_file,
        "-o", pangenome_manifest
    ]
    
    step = RunShellCommand(
        command=cmd,
        name="compile_pangenomes",
        validate_output_filepaths=[pangenome_manifest]
    ).run()
    step.check_status()
    logger.info(f"Pangenome manifest complete")
    logger.info("")
    
    # Only run protein clustering if proteins are provided
    has_proteins = any(data["proteins"] is not None for data in genome_data.values())
    
    if has_proteins:
        # ===============================
        # Step 4: Run protein clustering
        # ===============================
        logger.info("="*80)
        logger.info("Step 4: Protein clustering within pangenomes")
        logger.info("="*80)
        
        # Read manifest to show how many pangenomes will be processed
        df_pangenome_manifest = pd.read_csv(pangenome_manifest, sep="\t", header=None)
        n_pangenomes = df_pangenome_manifest.iloc[:, 1].nunique()
        logger.info(f"Processing {n_pangenomes} pangenomes")
        logger.info(f"(See {os.path.join(directories['project'], 'protein_clustering', 'log')} for detailed progress)")
        logger.info("")
        
        protein_clustering_dir = os.path.join(directories["project"], "protein_clustering")
        protein_clusters_file = os.path.join(protein_clustering_dir, "output", "proteins_to_orthologs.tsv.gz")
        
        cmd = [
            "pangenomium", "cluster-proteins-from-pangenomes",
            "-i", pangenome_manifest,
            "-o", protein_clustering_dir,
            "--n_threads_per_task", str(args.n_threads_mmseqs_per_task),
            "--n_concurrent_tasks", str(args.n_concurrent_mmseqs_tasks),
            "-a", args.algorithm,
            "-t", str(args.minimum_identity_threshold),
            "--minimum_coverage_threshold", str(args.minimum_coverage_threshold),
            "--separator", args.separator,
            "--protein_cluster_prefix", args.protein_cluster_prefix,
            "--cluster_label_mode", args.cluster_label_mode,
        ]
        
        if args.mmseqs2_options:
            cmd.extend(["--mmseqs2_options", args.mmseqs2_options])
        if args.no_singletons:
            cmd.append("--no_singletons")
        
        step = RunShellCommand(
            command=cmd,
            name="cluster_proteins",
            validate_output_filepaths=[protein_clusters_file]
        ).run()
        step.check_status()
        logger.info(f"Protein clustering complete")
        logger.info("")
    
    # ===============================
    # Step 5: Generate comprehensive outputs
    # ===============================
    logger.info("="*80)
    logger.info("Step 5: Generating comprehensive outputs")
    logger.info("="*80)
    
    generate_comprehensive_output(
        genome_data=genome_data,
        genome_clusters_file=genome_clusters_file,
        protein_clusters_file=protein_clusters_file if has_proteins else None,
        directories=directories,
        args=args,
    )
    
    # ===============================
    # Step 6: Cleanup temporary directories
    # ===============================
    if not args.keep_temporary:
        logger.info("="*80)
        logger.info("Step 6: Cleaning up temporary directories")
        logger.info("="*80)
        
        import shutil
        cleanup_dirs = []
        
        # Main project tmp and intermediate
        if os.path.exists(directories["tmp"]):
            cleanup_dirs.append(directories["tmp"])
        if os.path.exists(directories["intermediate"]):
            cleanup_dirs.append(directories["intermediate"])
        
        # Genome clustering tmp and intermediate
        genome_clustering_tmp = os.path.join(directories["project"], "genome_clustering", "tmp")
        genome_clustering_intermediate = os.path.join(directories["project"], "genome_clustering", "intermediate")
        if os.path.exists(genome_clustering_tmp):
            cleanup_dirs.append(genome_clustering_tmp)
        if os.path.exists(genome_clustering_intermediate):
            cleanup_dirs.append(genome_clustering_intermediate)
        
        # Protein clustering tmp and intermediate
        if has_proteins:
            protein_clustering_tmp = os.path.join(directories["project"], "protein_clustering", "tmp")
            protein_clustering_intermediate = os.path.join(directories["project"], "protein_clustering", "intermediate")
            if os.path.exists(protein_clustering_tmp):
                cleanup_dirs.append(protein_clustering_tmp)
            if os.path.exists(protein_clustering_intermediate):
                cleanup_dirs.append(protein_clustering_intermediate)
        
        # Remove directories
        for cleanup_dir in cleanup_dirs:
            logger.info(f"Removing: {cleanup_dir}")
            shutil.rmtree(cleanup_dir, ignore_errors=True)
        
        logger.info(f"Cleaned up {len(cleanup_dirs)} temporary directories")
        logger.info("")
    
    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    
    return 0


def generate_comprehensive_output(genome_data, genome_clusters_file, protein_clusters_file, directories, args):
    """Generate all VEBA-style comprehensive output files"""
    
    # Initialize tracking dictionaries
    genome_to_sample = {}
    genome_to_organism_type = {}
    genome_to_genomecluster = {}
    genome_to_number_of_contigs = defaultdict(int)
    genome_to_number_of_proteins = defaultdict(int)
    
    contig_to_genome = {}
    contig_to_sample = {}
    
    protein_to_genome = {}
    protein_to_sample = {}
    protein_to_sequence = {}
    protein_to_cds = {}
    protein_to_proteincluster = {}
    
    # =============================
    # Parse genome clusters
    # =============================
    logger.info("Loading genome clusters...")
    df_genome_clusters = pd.read_csv(genome_clusters_file, sep="\t", index_col=0, header=None)
    df_genome_clusters.columns = ["id_genome_cluster"]
    genome_to_genomecluster = df_genome_clusters["id_genome_cluster"].to_dict()
    
    # =============================
    # Parse genome FASTAs
    # =============================
    logger.info("Parsing genome FASTAs for contig information...")
    for id_genome, data in tqdm(genome_data.items(), desc="Processing genomes", unit="genome"):
        genome_to_sample[id_genome] = data["id_sample"]
        genome_to_organism_type[id_genome] = data["organism_type"]
        
        # Parse genome FASTA
        for id_contig, seq in pyfastx.Fasta(data["genome"], build_index=False):
            id_contig = id_contig.split()[0]  # Remove description
            contig_to_genome[id_contig] = id_genome
            contig_to_sample[id_contig] = data["id_sample"]
            genome_to_number_of_contigs[id_genome] += 1
    
    # =============================
    # Parse protein FASTAs
    # =============================
    if protein_clusters_file:
        logger.info("Parsing protein FASTAs...")
        for id_genome, data in tqdm(genome_data.items(), desc="Processing proteins", unit="genome"):
            if data["proteins"] is None:
                continue
            
            # Parse protein FASTA
            for id_protein, seq in pyfastx.Fasta(data["proteins"], build_index=False):
                id_protein = id_protein.split()[0]
                protein_to_genome[id_protein] = id_genome
                protein_to_sample[id_protein] = data["id_sample"]
                protein_to_sequence[id_protein] = seq
                genome_to_number_of_proteins[id_genome] += 1
            
            # Parse CDS FASTA if provided
            if data["cds"] is not None:
                for id_protein, seq in pyfastx.Fasta(data["cds"], build_index=False):
                    id_protein = id_protein.split()[0]
                    protein_to_cds[id_protein] = seq
        
        # Load protein clusters
        logger.info("Loading protein clusters...")
        df_protein_clusters = pd.read_csv(protein_clusters_file, sep="\t", index_col=0, header=None)
        df_protein_clusters.columns = ["id_protein_cluster"]
        protein_to_proteincluster = df_protein_clusters["id_protein_cluster"].to_dict()
    
    # Convert to Series
    genome_to_sample = pd.Series(genome_to_sample)
    genome_to_organism_type = pd.Series(genome_to_organism_type)
    genome_to_genomecluster = pd.Series(genome_to_genomecluster)
    genome_to_number_of_contigs = pd.Series(genome_to_number_of_contigs)
    genome_to_number_of_proteins = pd.Series(genome_to_number_of_proteins)
    
    contig_to_genome = pd.Series(contig_to_genome)
    contig_to_sample = pd.Series(contig_to_sample)
    
    if protein_clusters_file:
        protein_to_genome = pd.Series(protein_to_genome)
        protein_to_sample = pd.Series(protein_to_sample)
        protein_to_sequence = pd.Series(protein_to_sequence)
        protein_to_cds = pd.Series(protein_to_cds)
        protein_to_proteincluster = pd.Series(protein_to_proteincluster)
    
    # =============================
    # Create identifier mapping tables
    # =============================
    logger.info("Creating identifier mapping tables...")
    
    # Genomes
    df_genomes = pd.DataFrame(OrderedDict([
        ("organism_type", genome_to_organism_type),
        ("sample_of_origin", genome_to_sample),
        ("id_genome_cluster", genome_to_genomecluster),
        ("number_of_contigs", genome_to_number_of_contigs),
        ("number_of_proteins", genome_to_number_of_proteins),
    ])).sort_values(["sample_of_origin", "id_genome_cluster"])
    df_genomes.index.name = "id_genome"
    
    # Contigs
    df_contigs = pd.DataFrame(OrderedDict([
        ("organism_type", contig_to_genome.map(lambda x: genome_to_organism_type[x])),
        ("id_genome", contig_to_genome),
        ("sample_of_origin", contig_to_sample),
        ("id_genome_cluster", contig_to_genome.map(lambda x: genome_to_genomecluster.get(x))),
    ])).sort_values(["sample_of_origin", "id_genome", "id_genome_cluster"])
    df_contigs.index.name = "id_contig"
    
    # Proteins
    if protein_clusters_file:
        df_proteins = pd.DataFrame(OrderedDict([
            ("organism_type", protein_to_genome.map(lambda x: genome_to_organism_type[x])),
            ("id_genome", protein_to_genome),
            ("sample_of_origin", protein_to_sample),
            ("id_genome_cluster", protein_to_genome.map(lambda x: genome_to_genomecluster.get(x))),
            ("id_protein_cluster", protein_to_proteincluster),
        ])).sort_values(["sample_of_origin", "id_genome", "id_genome_cluster", "id_protein_cluster"])
        df_proteins.index.name = "id_protein"
    
    # =============================
    # Generate genome cluster summary
    # =============================
    logger.info("Generating genome cluster summary...")
    genomecluster_data = defaultdict(dict)
    for id_genomecluster, df in df_genomes.groupby("id_genome_cluster"):
        if pd.isnull(id_genomecluster):
            continue
        genomecluster_data[id_genomecluster]["number_of_components"] = df.shape[0]
        genomecluster_data[id_genomecluster]["components"] = set(df.index)
        genomecluster_data[id_genomecluster]["number_of_samples_of_origin"] = df["sample_of_origin"].nunique()
        genomecluster_data[id_genomecluster]["samples_of_origin"] = set(df["sample_of_origin"].unique())
    
    df_genomeclusters = pd.DataFrame(genomecluster_data).T.sort_index()
    df_genomeclusters = df_genomeclusters[["number_of_components", "number_of_samples_of_origin", "components", "samples_of_origin"]]
    df_genomeclusters.index.name = "id_genome_cluster"
    
    # =============================
    # Protein cluster analysis (if proteins provided)
    # =============================
    if protein_clusters_file:
        logger.info("Generating protein cluster summary...")
        
        # Protein cluster summary
        proteincluster_data = defaultdict(dict)
        for id_proteincluster, df in df_proteins.groupby("id_protein_cluster"):
            if pd.isnull(id_proteincluster):
                continue
            proteincluster_data[id_proteincluster]["number_of_components"] = df.shape[0]
            proteincluster_data[id_proteincluster]["components"] = set(df.index)
            proteincluster_data[id_proteincluster]["number_of_samples_of_origin"] = df["sample_of_origin"].nunique()
            proteincluster_data[id_proteincluster]["samples_of_origin"] = set(df["sample_of_origin"].unique())
        
        df_proteinclusters = pd.DataFrame(proteincluster_data).T.sort_index()
        df_proteinclusters = df_proteinclusters[["number_of_components", "number_of_samples_of_origin", "components", "samples_of_origin"]]
        df_proteinclusters.index.name = "id_protein_cluster"
        
        # Generate prevalence tables and analyze core/singletons
        logger.info("Generating prevalence tables and analyzing core pangenomes...")
        
        genome_to_number_of_singletons = []
        genome_to_ratio_of_singletons = []
        protein_to_number_of_genomes_detected = {}
        protein_to_ratio_of_genomes_detected = {}
        genomecluster_to_corepangenome = {}
        genomecluster_to_singletons = {}
        proteincluster_to_representative = {}
        
        for id_genomecluster, df in tqdm(df_proteins.groupby("id_genome_cluster"), 
                                        desc="Analyzing pangenomes", unit="pangenome"):
            if pd.isnull(id_genomecluster):
                continue
            
            # Prevalence matrix
            df_subset = df.reset_index()[["id_genome", "id_protein", "id_protein_cluster"]]
            df_prevalence = get_protein_cluster_prevalence(df_subset)
            
            # Write prevalence table
            prevalence_file = os.path.join(directories["output"], "pangenome_tables", f"{id_genomecluster}.tsv.gz")
            df_prevalence.to_csv(prevalence_file, sep="\t")
            
            # Detected/Not-detected analysis
            df_prevalence_binary = df_prevalence > 0
            number_of_genomes_detected = df_prevalence_binary.sum(axis=0)
            ratio_of_genomes_detected = df_prevalence_binary.mean(axis=0)
            protein_to_number_of_genomes_detected.update(number_of_genomes_detected.to_dict())
            protein_to_ratio_of_genomes_detected.update(ratio_of_genomes_detected.to_dict())
            
            # Core pangenome
            core_proteinclusters = ratio_of_genomes_detected[ratio_of_genomes_detected >= args.minimum_core_prevalence].index
            genomecluster_to_corepangenome[id_genomecluster] = set(core_proteinclusters)
            
            # Singletons
            singleton_proteinclusters = number_of_genomes_detected[number_of_genomes_detected == 1].index
            genomecluster_to_singletons[id_genomecluster] = set(singleton_proteinclusters)
            
            # Count singletons per genome (only for clusters with >1 genome)
            if df_prevalence.shape[0] > 1:
                singleton_prevalence = df_prevalence_binary.loc[:, singleton_proteinclusters]
                number_of_singletons = singleton_prevalence.sum(axis=1)
                ratio_of_singletons = number_of_singletons / df_prevalence_binary.sum(axis=1)
                genome_to_number_of_singletons.append(number_of_singletons)
                genome_to_ratio_of_singletons.append(ratio_of_singletons)
            
            # Get representatives (protein with highest connectivity in cluster)
            for id_proteincluster in df_subset["id_protein_cluster"].unique():
                if pd.isnull(id_proteincluster):
                    continue
                proteins_in_cluster = df_subset[df_subset["id_protein_cluster"] == id_proteincluster]["id_protein"]
                # For now, just use the first protein as representative
                # TODO: Could calculate connectivity if we had the graph
                proteincluster_to_representative[id_proteincluster] = proteins_in_cluster.iloc[0]
        
        # Add singleton statistics to genome table
        if genome_to_number_of_singletons:
            df_genomes["number_of_singleton_protein_clusters"] = pd.concat(genome_to_number_of_singletons).reindex(df_genomes.index).astype("Int64")
            df_genomes["ratio_of_protein_clusters_are_singletons"] = pd.concat(genome_to_ratio_of_singletons).reindex(df_genomes.index)
        
        # Add detection statistics to protein cluster table
        df_proteinclusters["number_of_genomes_detected"] = pd.Series(protein_to_number_of_genomes_detected)
        df_proteinclusters["ratio_of_genomes_detected"] = pd.Series(protein_to_ratio_of_genomes_detected)
        df_proteinclusters["core_pangenome"] = df_proteinclusters["ratio_of_genomes_detected"] >= args.minimum_core_prevalence
        df_proteinclusters["singleton"] = df_proteinclusters["number_of_genomes_detected"] == 1
        
        # Add core/singleton info to genome cluster table
        genomecluster_to_corepangenome = pd.Series(genomecluster_to_corepangenome)
        genomecluster_to_singletons = pd.Series(genomecluster_to_singletons)
        df_genomeclusters["number_of_proteins_in_core_pangenome"] = genomecluster_to_corepangenome.map(len)
        df_genomeclusters["core_pangenome"] = genomecluster_to_corepangenome
        df_genomeclusters["number_of_singleton_proteins"] = genomecluster_to_singletons.map(len)
        df_genomeclusters["singletons"] = genomecluster_to_singletons
        
        # Calculate average copies per genome
        number_of_gene_copies_per_sspc_per_genome = df_proteins.groupby(["id_genome", "id_protein_cluster"]).size()
        df_proteinclusters["average_number_of_copies_per_genome"] = number_of_gene_copies_per_sspc_per_genome.groupby(lambda x: x[1]).mean()
    
    # =============================
    # Calculate feature compression ratios
    # =============================
    if protein_clusters_file:
        logger.info("Calculating feature compression ratios...")
        fcr_data = defaultdict(dict)
        
        for organism_type, df in df_genomes.groupby("organism_type"):
            feature_to_cluster = df["id_genome_cluster"].dropna()
            fcr_data[organism_type]["number_of_genomes"] = feature_to_cluster.size
            fcr_data[organism_type]["number_of_genome_clusters"] = feature_to_cluster.nunique()
            fcr_data[organism_type]["genomic_fcr"] = 1 - (feature_to_cluster.nunique() / feature_to_cluster.size)
        
        for organism_type, df in df_proteins.groupby("organism_type"):
            feature_to_cluster = df["id_protein_cluster"].dropna()
            fcr_data[organism_type]["number_of_proteins"] = feature_to_cluster.size
            fcr_data[organism_type]["number_of_protein_clusters"] = feature_to_cluster.nunique()
            fcr_data[organism_type]["functional_fcr"] = 1 - (feature_to_cluster.nunique() / feature_to_cluster.size)
        
        df_fcr = pd.DataFrame(fcr_data).T.sort_index()
        df_fcr = df_fcr[["number_of_genomes", "number_of_genome_clusters", "genomic_fcr", 
                        "number_of_proteins", "number_of_protein_clusters", "functional_fcr"]]
        for field in ["number_of_genomes", "number_of_genome_clusters", "number_of_proteins", "number_of_protein_clusters"]:
            df_fcr[field] = df_fcr[field].astype(int)
        df_fcr.index.name = "organism_type"
    
    # =============================
    # Write representative sequences
    # =============================
    if protein_clusters_file and not args.no_representative_sequences:
        logger.info("Writing representative sequences...")
        with open_file_writer(os.path.join(directories["output"], "representatives", "representative_sequences.faa.gz")) as f:
            for id_proteincluster, id_representative in tqdm(proteincluster_to_representative.items(),
                                                            desc="Writing representatives", unit="cluster"):
                seq = protein_to_sequence[id_representative]
                print(f">{id_proteincluster} {id_representative}\n{seq}", file=f)
    
    # =============================
    # Write core pangenome sequences
    # =============================
    if protein_clusters_file and not args.no_core_sequences:
        logger.info("Writing core pangenome sequences...")
        for id_genomecluster, core_clusters in tqdm(genomecluster_to_corepangenome.items(),
                                                    desc="Writing core sequences", unit="pangenome"):
            # Protein sequences (gzipped)
            faa_file = os.path.join(directories["output"], "pangenome_core_sequences", f"{id_genomecluster}.faa.gz")
            with open_file_writer(faa_file) as f:
                for id_proteincluster in sorted(core_clusters):
                    if id_proteincluster not in proteincluster_to_representative:
                        continue
                    id_representative = proteincluster_to_representative[id_proteincluster]
                    seq = protein_to_sequence[id_representative]
                    print(f">{id_proteincluster} {id_representative}\n{seq}", file=f)
            
            # CDS sequences (if available, gzipped)
            if not protein_to_cds.empty:
                ffn_file = os.path.join(directories["output"], "pangenome_core_sequences", f"{id_genomecluster}.ffn.gz")
                with open_file_writer(ffn_file) as f:
                    for id_proteincluster in sorted(core_clusters):
                        if id_proteincluster not in proteincluster_to_representative:
                            continue
                        id_representative = proteincluster_to_representative[id_proteincluster]
                        if id_representative in protein_to_cds:
                            seq = protein_to_cds[id_representative]
                            print(f">{id_proteincluster} {id_representative}\n{seq}", file=f)
    
    # =============================
    # Write all output tables
    # =============================
    logger.info("Writing output tables...")
    
    # Identifier mappings
    df_genomes.to_csv(os.path.join(directories["output"], "identifier_mapping.genomes.tsv.gz"), sep="\t")
    df_contigs.to_csv(os.path.join(directories["output"], "identifier_mapping.contigs.tsv.gz"), sep="\t")
    if protein_clusters_file:
        df_proteins.to_csv(os.path.join(directories["output"], "identifier_mapping.proteins.tsv.gz"), sep="\t")
    
    # Cluster summaries
    df_genomeclusters.to_csv(os.path.join(directories["output"], "genome_clusters.tsv.gz"), sep="\t")
    if protein_clusters_file:
        df_proteinclusters.to_csv(os.path.join(directories["output"], "protein_clusters.tsv.gz"), sep="\t")
        df_fcr.to_csv(os.path.join(directories["output"], "feature_compression_ratios.tsv.gz"), sep="\t")
    
    # Simple mapping files
    df_genomes["id_genome_cluster"].to_frame().dropna(how="any", axis=0).to_csv(
        os.path.join(directories["output"], "genomes_to_pangenomes.tsv.gz"), sep="\t", header=None)
    df_contigs["id_genome"].to_frame().dropna(how="any", axis=0).to_csv(
        os.path.join(directories["output"], "contigs_to_genomes.tsv.gz"), sep="\t", header=None)
    df_contigs["id_genome_cluster"].to_frame().dropna(how="any", axis=0).to_csv(
        os.path.join(directories["output"], "contigs_to_pangenomes.tsv.gz"), sep="\t", header=None)
    if protein_clusters_file:
        df_proteins["id_protein_cluster"].to_frame().dropna(how="any", axis=0).to_csv(
            os.path.join(directories["output"], "proteins_to_orthologs.tsv.gz"), sep="\t", header=None)
    
    # Copy serialization files from clustering directories
    logger.info("Copying serialization files...")
    genome_clustering_dir = os.path.join(directories["project"], "genome_clustering")
    
    # For VEBA mode with organism types, use organism prefix
    if args.mode == "veba":
        # Group by organism type and copy files
        for organism_type in genome_to_organism_type.unique():
            org_prefix = organism_type.lower()
            
            # Genome serialization
            src_graph = os.path.join(genome_clustering_dir, "output", "serialization", "genome_clusters.graph.pkl.gz")
            src_dict = os.path.join(genome_clustering_dir, "output", "serialization", "genome_clusters.dict.pkl.gz")
            src_repr = os.path.join(genome_clustering_dir, "output", "representatives", "genome_representatives.tsv.gz")
            
            if os.path.exists(src_graph):
                dst_graph = os.path.join(directories["output"], "serialization", f"{org_prefix}.networkx_graph.pkl.gz")
                os.system(f"cp {src_graph} {dst_graph}")
            if os.path.exists(src_dict):
                dst_dict = os.path.join(directories["output"], "serialization", f"{org_prefix}.dict.pkl.gz")
                os.system(f"cp {src_dict} {dst_dict}")
            if os.path.exists(src_repr):
                dst_repr = os.path.join(directories["output"], "representatives", f"{org_prefix}.representatives.tsv.gz")
                os.system(f"cp {src_repr} {dst_repr}")
    else:
        # Batch mode: just copy without organism prefix
        src_graph = os.path.join(genome_clustering_dir, "output", "serialization", "genome_clusters.graph.pkl.gz")
        src_dict = os.path.join(genome_clustering_dir, "output", "serialization", "genome_clusters.dict.pkl.gz")
        src_repr = os.path.join(genome_clustering_dir, "output", "representatives", "genome_representatives.tsv.gz")
        
        if os.path.exists(src_graph):
            dst_graph = os.path.join(directories["output"], "serialization", "genome_clusters.graph.pkl.gz")
            os.system(f"cp {src_graph} {dst_graph}")
        if os.path.exists(src_dict):
            dst_dict = os.path.join(directories["output"], "serialization", "genome_clusters.dict.pkl.gz")
            os.system(f"cp {src_dict} {dst_dict}")
        if os.path.exists(src_repr):
            dst_repr = os.path.join(directories["output"], "representatives", "genome_representatives.tsv.gz")
            os.system(f"cp {src_repr} {dst_repr}")
    
    logger.info("Comprehensive output generation complete")
    logger.info(f"Output directory: {directories['output']}")
