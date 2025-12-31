"""Cluster proteins within pangenomes into pangenome-specific orthogroups"""
from __future__ import print_function, division
import sys
import os
import argparse
import gzip
import shutil
from multiprocessing import cpu_count
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
from tqdm import tqdm
from loguru import logger
from pyexeggutor import RunShellCommand, format_header, open_file_writer
from .. import __version__
from ..utils import setup_directories, setup_logger, print_header


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

def parse_input(args, directories):
    """Parse pangenome manifest and prepare protein files"""
    pangenome_to_proteins = defaultdict(dict)
    pangenome_to_protein_genome_map = defaultdict(dict)  # Track protein -> genome mapping
    
    if args.pangenome_manifest:
        # Direct pangenome manifest: [id_genome, id_pangenome, protein_filepath]
        df = pd.read_csv(args.pangenome_manifest, sep="\t", header=None)
        df.columns = ["id_genome", "id_pangenome", "protein_filepath"]
        
        for pangenome_id, group in df.groupby("id_pangenome"):
            # Create concatenated protein file for this pangenome (gzipped)
            protein_fasta = os.path.join(
                directories["intermediate"],
                f"{pangenome_id}_proteins.faa.gz"
            )
            
            protein_to_genome = {}
            
            with open_file_writer(protein_fasta) as f_out:
                for _, row in group.iterrows():
                    filepath = row["protein_filepath"]
                    genome_id = row["id_genome"]
                    
                    if filepath.endswith(".gz"):
                        opener = gzip.open
                        mode = "rt"
                    else:
                        opener = open
                        mode = "r"
                    
                    with opener(filepath, mode) as f_in:
                        for line in f_in:
                            if line.startswith(">"):
                                protein_id = line.strip()[1:].split()[0]
                                protein_to_genome[protein_id] = genome_id
                            f_out.write(line)
            
            pangenome_to_proteins[pangenome_id]["fasta"] = protein_fasta
            pangenome_to_proteins[pangenome_id]["genomes"] = group["id_genome"].tolist()
            pangenome_to_protein_genome_map[pangenome_id] = protein_to_genome
    
    elif args.genome_manifest and args.genome_clusters:
        # Create pangenome manifest from genome manifest + clusters
        df_genomes = pd.read_csv(args.genome_manifest, sep="\t", header=None)
        df_clusters = pd.read_csv(args.genome_clusters, sep="\t", header=None, names=["id_genome", "id_pangenome"])
        
        # Determine manifest format
        if df_genomes.shape[1] == 4:
            # Batch: [organism_type, id_genome, genome_filepath, protein_filepath]
            df_genomes.columns = ["organism_type", "id_genome", "genome_filepath", "protein_filepath"]
            df_genomes = df_genomes[["id_genome", "protein_filepath"]]
        elif df_genomes.shape[1] >= 5:
            # VEBA: [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, ...]
            df_genomes = df_genomes.iloc[:, [2, 4]]
            df_genomes.columns = ["id_genome", "protein_filepath"]
        
        # Merge with clusters
        df_merged = df_clusters.merge(df_genomes, on="id_genome", how="left")
        
        for pangenome_id, group in df_merged.groupby("id_pangenome"):
            protein_fasta = os.path.join(
                directories["intermediate"],
                f"{pangenome_id}_proteins.faa.gz"
            )
            
            protein_to_genome = {}
            
            with open_file_writer(protein_fasta) as f_out:
                for _, row in group.iterrows():
                    filepath = row["protein_filepath"]
                    genome_id = row["id_genome"]
                    
                    if pd.isna(filepath) or not os.path.exists(filepath):
                        continue
                    
                    if filepath.endswith(".gz"):
                        opener = gzip.open
                        mode = "rt"
                    else:
                        opener = open
                        mode = "r"
                    
                    with opener(filepath, mode) as f_in:
                        for line in f_in:
                            if line.startswith(">"):
                                protein_id = line.strip()[1:].split()[0]
                                protein_to_genome[protein_id] = genome_id
                            f_out.write(line)
            
            pangenome_to_proteins[pangenome_id]["fasta"] = protein_fasta
            pangenome_to_proteins[pangenome_id]["genomes"] = group["id_genome"].tolist()
            pangenome_to_protein_genome_map[pangenome_id] = protein_to_genome
    
    return pangenome_to_proteins, pangenome_to_protein_genome_map

def register_parser(subparsers):
    """Register the cluster-proteins-from-pangenomes subcommand parser"""
    parser = subparsers.add_parser(
        'cluster-proteins-from-pangenomes',
        help='Cluster proteins within pangenomes into pangenome-specific orthogroups',
        description='Cluster proteins within each pangenome separately using MMseqs2',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # I/O arguments
    parser_io = parser.add_argument_group('I/O arguments')
    parser_io.add_argument("-i", "--pangenome_manifest", type=str, help="Pangenome manifest: [id_genome, id_pangenome, protein_filepath]")
    parser_io.add_argument("-g", "--genome_manifest", type=str, help="Genome manifest (use with --genome_clusters)")
    parser_io.add_argument("-c", "--genome_clusters", type=str, help="Genome clusters file (use with --genome_manifest)")
    parser_io.add_argument("-o", "--output_directory", type=str, default="pangenome_protein_clustering_output", help="Output directory [Default: pangenome_protein_clustering_output]")
    
    # Utility arguments
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--n_threads_per_task", type=int, default=1, help="Threads per task [Default: 1]")
    parser_utility.add_argument("--n_concurrent_tasks", type=int, default=1, help="Number of tasks to process in parallel [Default: 1]")
    
    # MMseqs2 arguments
    parser_mmseqs = parser.add_argument_group('MMseqs2 arguments')
    parser_mmseqs.add_argument("-a", "--algorithm", type=str, default="mmseqs-cluster", choices=["mmseqs-cluster", "mmseqs-linclust"], help="Algorithm [Default: mmseqs-cluster]")
    parser_mmseqs.add_argument("-t", "--minimum_identity_threshold", type=float, default=50.0, help="Identity threshold [Default: 50.0]")
    parser_mmseqs.add_argument("--minimum_coverage_threshold", type=float, default=0.8, help="Coverage threshold [Default: 0.8]")
    parser_mmseqs.add_argument("--mmseqs2_options", type=str, default="", help="Additional MMseqs2 options")
    
    # Clustering arguments
    parser_clustering = parser.add_argument_group('Clustering arguments')
    parser_clustering.add_argument("--cluster_suffix", type=str, default="", help="Cluster suffix [Default: '']")
    parser_clustering.add_argument("--cluster_prefix_zfill", type=int, default=0, help="Prefix zfill [Default: 0]")
    parser_clustering.add_argument("--cluster_label_mode", type=str, default="md5", choices=["numeric", "random", "pseudo-random", "md5", "nodes"], help="Label mode [Default: md5]")
    parser_clustering.add_argument("--no_singletons", action="store_true", help="Exclude singletons")
    
    # Set the function to call
    parser.set_defaults(func=run)
    
    return parser

def process_pangenome(pangenome_id, protein_data, protein_to_genome, args, directories, algorithm):
    """Process a single pangenome - run MMseqs2, compile clusters, and create prevalence table"""
    
    protein_fasta = protein_data["fasta"]
    genomes = protein_data["genomes"]
    
    logger.info("="*80)
    logger.info(f"Clustering pangenome: {pangenome_id}")
    logger.info("="*80)
    
    # Create unique tmp directory for this pangenome to avoid conflicts in parallel execution
    pangenome_tmp = os.path.join(directories["tmp"], pangenome_id)
    os.makedirs(pangenome_tmp, exist_ok=True)
    
    # Run MMseqs2
    mmseqs_prefix = os.path.join(directories["intermediate"], f"{pangenome_id}_mmseqs2")
    edgelist_file = os.path.join(directories["intermediate"], f"{pangenome_id}_edgelist.tsv.gz")
    edgelist_temp = os.path.join(directories["intermediate"], f"{pangenome_id}_edgelist.tsv")
    
    cmd = [
        "mmseqs", algorithm,
        protein_fasta,
        mmseqs_prefix,
        pangenome_tmp,  # Use unique tmp directory
        "--threads", str(args.n_threads_per_task),
        "--min-seq-id", str(args.minimum_identity_threshold / 100),
        "-c", str(args.minimum_coverage_threshold),
        "--cov-mode", "1",
    ]
    
    if args.mmseqs2_options:
        cmd.append(args.mmseqs2_options)
    
    cmd.extend([
        "&&",
        "mv", f"{mmseqs_prefix}_cluster.tsv", edgelist_temp,
        "&&",
        "gzip", edgelist_temp,
        "&&",
        "rm", "-rf",
        f"{mmseqs_prefix}_all_seqs.fasta",
        f"{mmseqs_prefix}_rep_seq.fasta",
    ])
    
    step = RunShellCommand(
        command=cmd,
        name=f"mmseqs2_{pangenome_id}",
        validate_output_filepaths=[edgelist_file]
    ).run()
    step.check_status()
    
    # Compile clusters with pangenome ID as prefix
    # This creates cluster IDs like: PanG-abc123__PanOG-xyz789
    cluster_file = os.path.join(directories["output"], f"{pangenome_id}_protein_clusters.tsv.gz")
    
    cmd = [
        "edgelist-to-clusters.py",
        "-i", edgelist_file,
        "-o", cluster_file,
        "--cluster_prefix", f"{pangenome_id}__PanOG-",  # Use pangenome ID + PanOG- as prefix
        "--cluster_prefix_zfill", str(args.cluster_prefix_zfill),
        "--cluster_label_mode", args.cluster_label_mode,
    ]
    
    # Only add cluster_suffix if non-empty
    if args.cluster_suffix:
        cmd.extend(["--cluster_suffix", args.cluster_suffix])
    
    if args.no_singletons:
        cmd.append("--no_singletons")
    
    step = RunShellCommand(
        command=cmd,
        name=f"compile_{pangenome_id}",
        validate_output_filepaths=[cluster_file]
    ).run()
    step.check_status()
    
    # Create prevalence table (genomes × orthogroups)
    df_clusters = pd.read_csv(cluster_file, sep="\t", header=None, names=["id_protein", "id_protein_cluster"])
    df_clusters["id_genome"] = df_clusters["id_protein"].map(protein_to_genome)
    
    # Create prevalence matrix
    df_prevalence = get_protein_cluster_prevalence(df_clusters[["id_genome", "id_protein", "id_protein_cluster"]])
    
    # Save prevalence table
    prevalence_file = os.path.join(directories["output"], "pangenome_tables", f"{pangenome_id}.tsv.gz")
    df_prevalence.to_csv(prevalence_file, sep="\t")
    
    # Clean up pangenome-specific tmp directory
    if os.path.exists(pangenome_tmp):
        shutil.rmtree(pangenome_tmp, ignore_errors=True)
    
    logger.info("")
    return cluster_file


def run(args):
    """Execute cluster-proteins-from-pangenomes command"""
    
    # Validate
    if not args.pangenome_manifest and not (args.genome_manifest and args.genome_clusters):
        logger.error("Must provide either --pangenome_manifest OR both --genome_manifest and --genome_clusters")
        sys.exit(1)
    if args.pangenome_manifest and (args.genome_manifest or args.genome_clusters):
        logger.error("Cannot use --pangenome_manifest with --genome_manifest/--genome_clusters")
        sys.exit(1)
    
    if args.n_threads_per_task == -1:
        args.n_threads_per_task = cpu_count()
    if args.n_concurrent_tasks == -1:
        args.n_concurrent_tasks = cpu_count()
    
    # Setup directories with command-specific subdirectory
    directories = setup_directories(args.output_directory, subdirectory="pangenome_protein_clustering")
    
    # Create pangenome_tables subdirectory for prevalence matrices
    os.makedirs(os.path.join(directories["output"], "pangenome_tables"), exist_ok=True)
    
    # Setup logger
    setup_logger(directories["log"], "cluster_proteins_from_pangenomes.log")
    
    # Print info
    logger.info("="*80)
    logger.info("pangenomium cluster-proteins-from-pangenomes")
    logger.info("="*80)
    print_header(
        version=__version__, 
        n_jobs=f"{args.n_concurrent_tasks} concurrent tasks × {args.n_threads_per_task} threads per task"
    )
    
    # Parse input
    logger.info("-"*80)
    logger.info("Parsing input")
    logger.info("-"*80)
    pangenome_to_proteins, pangenome_to_protein_genome_map = parse_input(args, directories)
    logger.info(f"Pangenomes: {len(pangenome_to_proteins)}")
    logger.info("")
    
    # Process each pangenome
    all_cluster_files = []
    algorithm = "easy-cluster" if args.algorithm == "mmseqs-cluster" else "easy-linclust"
    
    # Use parallel execution if n_concurrent_tasks > 1
    if args.n_concurrent_tasks > 1:
        logger.info(f"Processing {len(pangenome_to_proteins)} pangenomes with {args.n_concurrent_tasks} concurrent tasks")
        logger.info("")
        
        with ThreadPoolExecutor(max_workers=args.n_concurrent_tasks) as executor:
            # Submit all jobs
            future_to_pangenome = {
                executor.submit(
                    process_pangenome, 
                    pangenome_id, 
                    protein_data, 
                    pangenome_to_protein_genome_map[pangenome_id],
                    args, 
                    directories, 
                    algorithm
                ): pangenome_id
                for pangenome_id, protein_data in pangenome_to_proteins.items()
            }
            
            # Process completed jobs with progress bar
            for future in tqdm(as_completed(future_to_pangenome), 
                             total=len(pangenome_to_proteins),
                             desc="Processing pangenomes", 
                             unit=" pangenomes"):
                pangenome_id = future_to_pangenome[future]
                try:
                    cluster_file = future.result()
                    all_cluster_files.append(cluster_file)
                except Exception as exc:
                    logger.error(f"Pangenome {pangenome_id} generated an exception: {exc}")
                    raise
    else:
        # Sequential execution (original behavior)
        for pangenome_id, protein_data in tqdm(pangenome_to_proteins.items(), desc="Processing pangenomes", unit=" pangenomes"):
            cluster_file = process_pangenome(
                pangenome_id, 
                protein_data, 
                pangenome_to_protein_genome_map[pangenome_id],
                args, 
                directories, 
                algorithm
            )
            all_cluster_files.append(cluster_file)
    
    # Concatenate all cluster files
    if len(all_cluster_files) > 0:
        logger.info("="*80)
        logger.info("Concatenating cluster files")
        logger.info("="*80)
        final_output = os.path.join(directories["output"], "proteins_to_orthologs.tsv.gz")
        
        cmd = ["cat"] + all_cluster_files + [">", final_output]
        
        step = RunShellCommand(
            command=cmd,
            name="concatenate",
            validate_output_filepaths=[final_output]
        ).run()
        step.check_status()
        
        # Clean up per-pangenome cluster files (info is in proteins_to_orthologs.tsv.gz)
        logger.info("Cleaning up per-pangenome cluster files...")
        for cluster_file in all_cluster_files:
            if os.path.exists(cluster_file):
                os.remove(cluster_file)
        
        logger.info("")
    
    logger.info("="*80)
    logger.info("Complete")
    logger.info("="*80)
    return 0
