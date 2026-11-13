# Configuration of the subset test (sourced by run_subset.sh). Every path is read only, except OUT.
# Defaults: the lab cluster (hn01, PBS), where the stored inputs and results of the paper live. The folders are
# readable by members of the group bliu.

## where the test writes (created; must not exist or be empty)
OUT=${OUT:-$HOME/bohr_subset_test}

## the subset: tissues and one chromosome (small tissues and chr21 keep the run short)
TISSUES=${TISSUES:-"Kidney_Cortex Brain_Substantia_nigra"}
CHR=${CHR:-chr21}
NCORES=${NCORES:-8}

## software
RSCRIPT=${RSCRIPT:-/home/e1327989/.local/share/mamba/envs/gtexfig/bin/Rscript}   # R 4.5 with data.table, VariantAnnotation, Matrix
BOHR_SRC=${BOHR_SRC:-}            # the BOHR package source (repository root or .tar.gz); empty = the copy that contains this paper/ folder

## stored inputs and results of the paper (GTEx v11; individual-level, dbGaP-protected: never copy them elsewhere)
GTEX=${GTEX:-/home/e1327989/project/GTEx}                                    # v11/ (proActiv outputs), metadata/, peer/v11/
VCF=${VCF:-/home/e1327989/project/GTEx/Rawdata/v9_vcf_MAF_0.05}              # chr<N>.vcf.gz + .tbi (GTEx v9 genotypes)
RUN_ROOT=${RUN_ROOT:-/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2} # work/final (tissue inputs), perm_maps, results/

## PBS (run_subset.sh submits itself when it is not already inside a job; SUBMIT=0 runs in the current shell)
SUBMIT=${SUBMIT:-1}
PBS_QUEUE=${PBS_QUEUE:-long}            # lab cluster: normal allows at most 24 GB, long 48 GB
PBS_RESOURCES=${PBS_RESOURCES:-select=1:ncpus=8:mem=48gb}
PBS_WALLTIME=${PBS_WALLTIME:-06:00:00}
PBS_PROJECT=${PBS_PROJECT:-}       # e.g. 11003054 on NSCC
