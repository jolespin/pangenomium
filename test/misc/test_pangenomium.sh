#!/bin/bash
# Comprehensive testing script for pangenomium v2025.1.15
# Tests all modules with all input formats

set -e  # Exit on error

echo "=================================================="
echo "Pangenomium Comprehensive Testing"
echo "=================================================="
echo ""

# Create directories
mkdir -p test_inputs test_outputs

# ============================================
# STEP 1: Generate different input formats
# ============================================
echo "Step 1: Generating test input files..."
echo ""

# 1a. Simple genome list (for cluster-genomes)
cut -f4 genomes_table.tsv > test_inputs/genome_list.txt
echo "✓ Created: test_inputs/genome_list.txt (simple genome list)"

# 1b. Simple protein list (for cluster-proteins)
cut -f5 genomes_table.tsv > test_inputs/protein_list.txt
echo "✓ Created: test_inputs/protein_list.txt (simple protein list)"

# 1c. Batch manifest (4 columns: organism_type, id_genome, genome_filepath, protein_filepath)
cut -f1,3,4,5 genomes_table.tsv > test_inputs/batch_manifest.tsv
echo "✓ Created: test_inputs/batch_manifest.tsv (batch format)"

# 1d. VEBA manifest (original - 7 columns)
cp genomes_table.tsv test_inputs/veba_manifest.tsv
echo "✓ Created: test_inputs/veba_manifest.tsv (VEBA format)"

echo ""
echo "Input files created in test_inputs/"
echo ""

# ============================================
# STEP 2: Test cluster-genomes with all modes
# ============================================
echo "=================================================="
echo "Step 2: Testing cluster-genomes"
echo "=================================================="
echo ""

# 2a. Simple list mode
echo "2a. cluster-genomes with simple list..."
pangenomium cluster-genomes \
  -i test_inputs/genome_list.txt \
  -o test_outputs/01_cluster_genomes_simple \
  --n_threads 4 \
  --ani_threshold 95.0
echo "✓ Completed: test_outputs/01_cluster_genomes_simple"
echo ""

# 2b. Batch mode
echo "2b. cluster-genomes with batch manifest..."
pangenomium cluster-genomes \
  -i test_inputs/batch_manifest.tsv \
  -o test_outputs/02_cluster_genomes_batch \
  --n_threads 4 \
  --ani_threshold 95.0
echo "✓ Completed: test_outputs/02_cluster_genomes_batch"
echo ""

# 2c. VEBA mode
echo "2c. cluster-genomes with VEBA manifest..."
pangenomium cluster-genomes \
  -i test_inputs/veba_manifest.tsv \
  -o test_outputs/03_cluster_genomes_veba \
  --n_threads 4 \
  --ani_threshold 95.0
echo "✓ Completed: test_outputs/03_cluster_genomes_veba"
echo ""

# 2d. Custom genome IDs (different from filenames)
echo "2d. cluster-genomes with custom genome IDs..."

# Create a 3-column manifest where genome IDs don't match filenames
# Format: [organism_type, id_genome, genome_filepath]
awk -F'\t' 'BEGIN{OFS="\t"; n=1} {
    printf "%s\tGENOME_%03d\t%s\n", $1, n, $4;
    n++
}' genomes_table.tsv > test_inputs/custom_id_manifest.tsv

pangenomium cluster-genomes \
  -i test_inputs/custom_id_manifest.tsv \
  -o test_outputs/02d_cluster_genomes_custom_ids \
  --n_threads 4 \
  --ani_threshold 95.0
echo "✓ Completed: test_outputs/02d_cluster_genomes_custom_ids"
echo ""

# ============================================
# STEP 3: Test cluster-proteins
# ============================================
echo "=================================================="
echo "Step 3: Testing cluster-proteins"
echo "=================================================="
echo ""

# 3a. Concatenate all proteins into single file
echo "3a. Creating concatenated protein file..."
cat $(cut -f5 genomes_table.tsv) > test_inputs/all_proteins.faa
echo "✓ Created: test_inputs/all_proteins.faa"
echo ""

# 3b. Run cluster-proteins
echo "3b. cluster-proteins..."
pangenomium cluster-proteins \
  -i test_inputs/all_proteins.faa \
  -o test_outputs/04_cluster_proteins \
  --n_threads 4 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/04_cluster_proteins"
echo ""

# ============================================
# STEP 4: Test cluster-proteins-from-pangenomes
# ============================================
echo "=================================================="
echo "Step 4: Testing cluster-proteins-from-pangenomes"
echo "=================================================="
echo ""

# 4a. With genome_manifest + genome_clusters (from batch)
echo "4a. cluster-proteins-from-pangenomes with genome_manifest + genome_clusters (batch)..."
pangenomium cluster-proteins-from-pangenomes \
  -g test_inputs/batch_manifest.tsv \
  -c test_outputs/02_cluster_genomes_batch/output/genomes_to_pangenomes.tsv.gz \
  -o test_outputs/05_cluster_proteins_from_pangenomes_batch \
  --n_threads_per_task 2 \
  --n_concurrent_tasks 4 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/05_cluster_proteins_from_pangenomes_batch"
echo ""

# 4b. With genome_manifest + genome_clusters (from VEBA)
echo "4b. cluster-proteins-from-pangenomes with genome_manifest + genome_clusters (VEBA)..."
pangenomium cluster-proteins-from-pangenomes \
  -g test_inputs/veba_manifest.tsv \
  -c test_outputs/03_cluster_genomes_veba/output/genomes_to_pangenomes.tsv.gz \
  -o test_outputs/06_cluster_proteins_from_pangenomes_veba \
  --n_threads_per_task 2 \
  --n_concurrent_tasks 4 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/06_cluster_proteins_from_pangenomes_veba"
echo ""

# 4c. Create pangenome_manifest and test with -i (from VEBA workflow)
echo "4c. Creating pangenome_manifest..."
# Format: [id_genome, id_pangenome, protein_filepath]
# Join genome clusters with genome manifest to get protein paths
join -t $'\t' \
  <(gunzip -c test_outputs/03_cluster_genomes_veba/output/genomes_to_pangenomes.tsv.gz | sort -k1,1) \
  <(awk -F'\t' '{print $3"\t"$5}' test_inputs/veba_manifest.tsv | sort -k1,1) | \
  awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3}' > test_inputs/pangenome_manifest.tsv

echo "✓ Created: test_inputs/pangenome_manifest.tsv"
echo ""

echo "4d. cluster-proteins-from-pangenomes with pangenome_manifest..."
pangenomium cluster-proteins-from-pangenomes \
  -i test_inputs/pangenome_manifest.tsv \
  -o test_outputs/07_cluster_proteins_from_pangenomes_direct \
  --n_threads_per_task 2 \
  --n_concurrent_tasks 4 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/07_cluster_proteins_from_pangenomes_direct"
echo ""

# ============================================
# STEP 5: Test end-to-end workflows
# ============================================
echo "=================================================="
echo "Step 5: Testing end-to-end workflows"
echo "=================================================="
echo ""

# 5a. end-to-end with batch mode
echo "5a. end-to-end with batch mode..."
pangenomium end-to-end \
  -i test_inputs/batch_manifest.tsv \
  --mode batch \
  -o test_outputs/08_end_to_end_batch \
  --n_threads_skani 4 \
  --n_threads_mmseqs_per_task 2 \
  --n_concurrent_mmseqs_tasks 4 \
  --ani_threshold 95.0 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/08_end_to_end_batch"
echo ""

# 5b. end-to-end with VEBA mode
echo "5b. end-to-end with VEBA mode..."
pangenomium end-to-end \
  -i test_inputs/veba_manifest.tsv \
  --mode veba \
  -o test_outputs/09_end_to_end_veba \
  --n_threads_skani 4 \
  --n_threads_mmseqs_per_task 2 \
  --n_concurrent_mmseqs_tasks 4 \
  --ani_threshold 95.0 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/09_end_to_end_veba"
echo ""

# ============================================
# STEP 6: Test manual end-to-end (step-by-step)
# ============================================
echo "=================================================="
echo "Step 6: Testing manual end-to-end workflow"
echo "=================================================="
echo ""

# 6a. cluster-genomes
echo "6a. Manual step 1: cluster-genomes..."
pangenomium cluster-genomes \
  -i test_inputs/veba_manifest.tsv \
  -o test_outputs/10_manual_step1_genomes \
  --n_threads 4 \
  --ani_threshold 95.0
echo "✓ Completed: test_outputs/10_manual_step1_genomes"
echo ""

# 6b. cluster-proteins-from-pangenomes
echo "6b. Manual step 2: cluster-proteins-from-pangenomes..."
pangenomium cluster-proteins-from-pangenomes \
  -g test_inputs/veba_manifest.tsv \
  -c test_outputs/10_manual_step1_genomes/output/genomes_to_pangenomes.tsv.gz \
  -o test_outputs/11_manual_step2_proteins \
  --n_threads_per_task 2 \
  --n_concurrent_tasks 4 \
  --minimum_identity_threshold 50.0
echo "✓ Completed: test_outputs/11_manual_step2_proteins"
echo ""

echo "=================================================="
echo "All tests completed successfully!"
echo "=================================================="
echo ""

# ============================================
# STEP 7: Verification checks
# ============================================
echo "=================================================="
echo "Step 7: Verification"
echo "=================================================="
echo ""

echo "Checking outputs..."
echo ""

# Check genome clustering outputs
echo "Genome clustering outputs:"
for dir in test_outputs/01_cluster_genomes_simple \
           test_outputs/02_cluster_genomes_batch \
           test_outputs/03_cluster_genomes_veba; do
  if [ -f "$dir/output/genomes_to_pangenomes.tsv.gz" ]; then
    n_clusters=$(gunzip -c "$dir/output/genomes_to_pangenomes.tsv.gz" | cut -f2 | sort -u | wc -l)
    echo "  ✓ $dir: $n_clusters pangenomes"
  fi
done
echo ""

# Check custom genome ID output
echo "Custom genome ID clustering:"
dir="test_outputs/02d_cluster_genomes_custom_ids"
if [ -f "$dir/output/genomes_to_pangenomes.tsv.gz" ]; then
  n_clusters=$(gunzip -c "$dir/output/genomes_to_pangenomes.tsv.gz" | cut -f2 | sort -u | wc -l)
  n_custom=$(gunzip -c "$dir/output/genomes_to_pangenomes.tsv.gz" | cut -f1 | grep -c "GENOME_")
  n_total=$(gunzip -c "$dir/output/genomes_to_pangenomes.tsv.gz" | wc -l)
  echo "  ✓ $dir: $n_clusters pangenomes, $n_custom/$n_total genomes have custom IDs"
  if [ "$n_custom" -eq "$n_total" ]; then
    echo "  ✓ All genome IDs are custom (symlink approach works!)"
  else
    echo "  ✗ FAIL: Some genome IDs are NOT custom"
    exit 1
  fi
fi
echo ""

# Check protein clustering from pangenomes
echo "Protein clustering outputs:"
for dir in test_outputs/05_cluster_proteins_from_pangenomes_batch \
           test_outputs/06_cluster_proteins_from_pangenomes_veba \
           test_outputs/07_cluster_proteins_from_pangenomes_direct; do
  if [ -f "$dir/protein_clustering/output/proteins_to_orthologs.tsv.gz" ]; then
    n_orthologs=$(gunzip -c "$dir/protein_clustering/output/proteins_to_orthologs.tsv.gz" | cut -f2 | sort -u | wc -l)
    n_tables=$(ls "$dir/protein_clustering/output/pangenome_tables/"*.tsv.gz 2>/dev/null | wc -l)
    echo "  ✓ $dir: $n_orthologs orthogroups, $n_tables pangenome tables"
  fi
done
echo ""

# Check end-to-end outputs
echo "End-to-end outputs:"
for dir in test_outputs/08_end_to_end_batch \
           test_outputs/09_end_to_end_veba; do
  if [ -f "$dir/output/protein_clusters.tsv.gz" ]; then
    n_genomes=$(gunzip -c "$dir/output/genome_clusters.tsv.gz" | tail -n +2 | wc -l)
    n_proteins=$(gunzip -c "$dir/output/protein_clusters.tsv.gz" | tail -n +2 | wc -l)
    n_tables=$(ls "$dir/output/pangenome_tables/"*.tsv.gz 2>/dev/null | wc -l)
    echo "  ✓ $dir: $n_genomes genomes, $n_proteins proteins, $n_tables pangenome tables"
  fi
done
echo ""

# Check manual workflow
echo "Manual workflow outputs:"
if [ -f "test_outputs/10_manual_step1_genomes/output/genomes_to_pangenomes.tsv.gz" ] && \
   [ -f "test_outputs/11_manual_step2_proteins/output/proteins_to_orthologs.tsv.gz" ]; then
  n_pangenomes=$(gunzip -c test_outputs/10_manual_step1_genomes/output/genomes_to_pangenomes.tsv.gz | cut -f2 | sort -u | wc -l)
  n_tables=$(ls test_outputs/11_manual_step2_proteins/output/pangenome_tables/*.tsv.gz 2>/dev/null | wc -l)
  echo "  ✓ Manual workflow: $n_pangenomes pangenomes, $n_tables pangenome tables"
fi
echo ""

# Compare manual vs end-to-end
echo "=================================================="
echo "Comparing manual workflow vs end-to-end"
echo "=================================================="
echo ""

echo "Genome clusters comparison:"
echo "  Manual:      $(gunzip -c test_outputs/10_manual_step1_genomes/output/genomes_to_pangenomes.tsv.gz | wc -l) lines"
echo "  End-to-end:  $(gunzip -c test_outputs/09_end_to_end_veba/genome_clustering/output/genomes_to_pangenomes.tsv.gz | wc -l) lines"
echo ""

echo "Protein clusters comparison:"
echo "  Manual:      $(gunzip -c test_outputs/11_manual_step2_proteins/output/proteins_to_orthologs.tsv.gz | wc -l) lines"
echo "  End-to-end:  $(gunzip -c test_outputs/09_end_to_end_veba/protein_clustering/output/proteins_to_orthologs.tsv.gz | wc -l) lines"
echo ""

echo "Pangenome tables comparison:"
echo "  Manual:      $(ls test_outputs/11_manual_step2_proteins/output/pangenome_tables/*.tsv.gz 2>/dev/null | wc -l) tables"
echo "  End-to-end:  $(ls test_outputs/09_end_to_end_veba/output/pangenome_tables/*.tsv.gz | wc -l) tables"
echo ""

echo "=================================================="
echo "Testing complete!"
echo "=================================================="
