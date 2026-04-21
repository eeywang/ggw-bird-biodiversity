#!/bin/bash
# =============================================================================
# submit_all.sh — Example SLURM submission scripts for the GGW bird pipeline
#
# Usage: Edit COUNTRY and N_SPECIES, then submit with:
#   sbatch submit_all.sh           (won't work — use individual sections below)
#
# Or submit each stage individually by extracting the sbatch blocks.
# =============================================================================

# ============================
# EDIT THESE FOR EACH COUNTRY
# ============================
COUNTRY="Ethiopia"           # Nigeria | Senegal | Ethiopia
N_SPECIES=477                # Number of species with ≥30 records
SCRATCH_DIR="/mnt/scratch/eeywang/ggw_birds/${COUNTRY}"

echo "=== GGW Bird Pipeline: ${COUNTRY} ==="
echo "Species: ${N_SPECIES}"
echo "Directory: ${SCRATCH_DIR}"

# =============================================================================
# STAGE 1: SDM Training + Hindcast (SLURM array)
# =============================================================================
# Save this block as: biomod2_job.sh
cat << 'EOF' > ${SCRATCH_DIR}/scripts/biomod2_job.sh
#!/bin/bash
#SBATCH --job-name=ggw_biomod2
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-477%50
#SBATCH --partition=nodes

module load miniforge
conda activate r-biomod2

cd /mnt/scratch/eeywang/ggw_birds/Ethiopia
mkdir -p logs

Rscript scripts/02_biomod2_modeling.R $SLURM_ARRAY_TASK_ID Ethiopia
EOF

# =============================================================================
# STAGE 2: Scenario Projection S1/S2 (SLURM array)
# =============================================================================
cat << 'EOF' > ${SCRATCH_DIR}/scripts/scenario_job.sh
#!/bin/bash
#SBATCH --job-name=ggw_scenario
#SBATCH --output=logs/scenario_%a.log
#SBATCH --array=1-470%50
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=nodes

module load miniforge
conda activate r-biomod2

cd /mnt/scratch/eeywang/ggw_birds/Ethiopia
Rscript scripts/01_scenario_projection.R $SLURM_ARRAY_TASK_ID Ethiopia
EOF

# =============================================================================
# STAGE 3: Richness Calculation (single job, ≥50 threshold)
# =============================================================================
cat << 'EOF' > ${SCRATCH_DIR}/scripts/richness_job.sh
#!/bin/bash
#SBATCH --job-name=ggw_richness
#SBATCH --output=logs/richness_%j.out
#SBATCH --error=logs/richness_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

module load miniforge
conda activate r-biomod2

cd /mnt/scratch/eeywang/ggw_birds/Ethiopia
mkdir -p logs
Rscript scripts/03_calculate_richness.R Ethiopia
EOF

# =============================================================================
# STAGE 4: Scenario Richness + Decomposition (single job, ≥50 threshold)
# =============================================================================
cat << 'EOF' > ${SCRATCH_DIR}/scripts/scenario_richness_job.sh
#!/bin/bash
#SBATCH --job-name=ggw_scen_rich
#SBATCH --output=logs/scenario_richness_%j.out
#SBATCH --error=logs/scenario_richness_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

module load miniforge
conda activate r-biomod2

cd /mnt/scratch/eeywang/ggw_birds/Ethiopia
mkdir -p logs
Rscript scripts/02_scenario_richness.R Ethiopia
EOF

# =============================================================================
# STAGE 5: Guild Richness (single job)
# =============================================================================
cat << 'EOF' > ${SCRATCH_DIR}/scripts/guild_job.sh
#!/bin/bash
#SBATCH --job-name=ggw_guild
#SBATCH --output=logs/guild_%j.out
#SBATCH --error=logs/guild_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=nodes

module load miniforge
conda activate r-biomod2

cd /mnt/scratch/eeywang/ggw_birds/Ethiopia
mkdir -p logs
Rscript scripts/03_guild_richness.R Ethiopia
EOF

echo ""
echo "Job scripts created in ${SCRATCH_DIR}/scripts/"
echo ""
echo "Submit in order:"
echo "  cd ${SCRATCH_DIR}/scripts"
echo "  sbatch biomod2_job.sh          # wait until all complete"
echo "  sbatch scenario_job.sh         # wait until all complete"
echo "  sbatch richness_job.sh"
echo "  sbatch scenario_richness_job.sh"
echo "  sbatch guild_job.sh"
