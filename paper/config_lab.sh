# Configuration of run_paper.sh (authored for the repository; sourced by run_paper.sh). Every path is read only, except
# OUT. Defaults: the lab cluster (hn01, PBS), where the stored inputs and results of the paper live, readable by members of
# the group bliu. A value already set in the environment wins over the default here.

## where the run writes (created; use a new or empty folder; delete it when done, it holds per-sample values)
OUT=${OUT:-$HOME/bohr_paper_run}

## what runs
STEPS=${STEPS:-"00 01 02 03 04 05 06 07"}   # add P for the package route (paper/with_package/)
ONLY=${ONLY:-}                              # optional list of sub-steps, e.g. "01.2 03.1"; empty = every sub-step of STEPS
KEEP_GOING=${KEEP_GOING:-0}                 # 1 = continue after a failed check
RESUME=${RESUME:-1}                         # 1 = skip sub-steps that passed in an earlier run into the same OUT

## the subset of the re-runs
TISSUES=${TISSUES:-"Kidney_Cortex Brain_Substantia_nigra"}   # 103 and 183 donors
CHR=${CHR:-chr21}
TRAITS=${TRAITS:-"SEX AGE"}                 # donor characteristics re-run in each tissue (step 04)
SIM_N=${SIM_N:-50}                          # simulated sample size re-run (step 02)
SIM_CHUNK=${SIM_CHUNK:-"1 4"}               # chunk k of n of its 2,000 simulated promoters fitted again (step 02)
NCORES=${NCORES:-8}

## software
RSCRIPT=${RSCRIPT:-/home/e1327989/.local/share/mamba/envs/gtexfig/bin/Rscript}   # R 4.5: data.table, VariantAnnotation, Matrix, coloc 5, gamlss
PYTHON=${PYTHON:-/home/e1327989/.local/share/mamba/envs/gtexfig/bin/python}      # Python 3 with pandas and pyarrow
SIM_RLIB=${SIM_RLIB:-/projects/group/bliu/e1327989/GTEx/new_simulation/newmethod_bohr2/m1_trunc/Rlib}   # truncnorm, as in the simulation runs
BOHR_SRC=${BOHR_SRC:-}                      # the BOHR package source; empty = the repository that holds this paper/ folder

## stored inputs and results of the paper (GTEx v11 individual-level data are dbGaP-protected: never copy them elsewhere)
GTEX=${GTEX:-/home/e1327989/project/GTEx}                                     # v11/ (proActiv outputs), metadata/, fit/, version_comparison/
VCF=${VCF:-/home/e1327989/project/GTEx/Rawdata/v9_vcf_MAF_0.05}               # chr<N>.vcf.gz + .tbi (GTEx v9 genotypes)
RUN_ROOT=${RUN_ROOT:-/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2}  # cis mapping, units, conditional analysis, eQTL coloc
SIM_ROOT=${SIM_ROOT:-/projects/group/bliu/e1327989/GTEx/new_simulation}        # simulations and benchmark (newmethod_bohr2/, figures_bohr2/)
TRG_ROOT=${TRG_ROOT:-/projects/group/bliu/e1327989/GTEx/tissue_regression}     # donor-trait runs
GWAS_ROOT=${GWAS_ROOT:-/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3}   # GWAS coloc, MRAS
F5S_V6=${F5S_V6:-/home/e1327989/f5_stage_20261106_v6}                         # major/minor stage (in/, out/)
F5S_V6C=${F5S_V6C:-/home/e1327989/f5_stage_20261107_v6c}                      # the same, three-promoter rule
F5_ENCODE=${F5_ENCODE:-/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v3}   # outputs of the ENCODE scripts
AG_ROOT=${AG_ROOT:-/projects/group/bliu/e1327989}                               # alphagenome_results/ (stored predictions) and GTEx/alphagenome/
SUPP=${SUPP:-/projects/group/bliu/e1327989/GTEx/manuscript/BOHR_v4/tables/supplementary}   # the published Supplementary Tables, one TSV per sheet

## PBS (run_paper.sh submits itself when it is not already inside a job; SUBMIT=0 runs in the current shell)
SUBMIT=${SUBMIT:-1}
PBS_QUEUE=${PBS_QUEUE:-long}                # lab cluster: normal allows at most 24 GB, long 48 GB
PBS_RESOURCES=${PBS_RESOURCES:-select=1:ncpus=8:mem=48gb}
PBS_WALLTIME=${PBS_WALLTIME:-12:00:00}
PBS_PROJECT=${PBS_PROJECT:-}                # e.g. 11003054 on NSCC
