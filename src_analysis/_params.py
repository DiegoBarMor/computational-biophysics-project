from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////
####################################### PARAMETERS

DIR_PROJECT = Path.cwd().parent
RUN_DETAILED_ANALYSIS = "mt2_rep0"
RUN_PREFFIXES = ["mt1", "mt2", "wt1", "wt2"]
RUNS = ["mt1_rep0", "mt1_rep1", "mt2_rep0", "mt2_rep1", "wt1_rep0", "wt1_rep1", "wt2_rep0", "wt2_rep1"]

CLUSTERING_LINK_METHOD = "ward"
COLOR_MAP_CLUSTERS = "prism"

N_RESIDS_SUBUNIT = 306 # hardcoded for systems 8DFN and 7SI9

####################################### DIRECTORIES

DIR_DA = DIR_PROJECT / "data_analysis"
DIR_DEXTRA = DIR_PROJECT / "data_extra"

DIR_DA_COORDS = DIR_DA / "_coords"
DIR_DA_TRAJECTORIES = DIR_DA / "_trajectories"
DIR_DA_BSE = DIR_DA / "bse"
DIR_DA_CLUSTERING = DIR_DA / "clustering"
DIR_DA_PCA = DIR_DA / "pca"
DIR_DA_PYINTERAPH = DIR_DA / "pyinteraph"
DIR_DA_RAMA = DIR_DA / "rama"
DIR_DA_RGYR = DIR_DA / "rgyr"
DIR_DA_RMSD = DIR_DA / "rmsd"
DIR_DA_RMSF = DIR_DA / "rmsf"
DIR_DA_SASA = DIR_DA / "sasa"
DIR_DA_VMD = DIR_DA / "vmd"

####################################### COLORS

RED   = 1, 0, 0, 1
GREEN = 0, 1, 0, 1
BLUE  = 0, 0, 1, 1

HALF_RED   = 1, 0, 0, 0.5
HALF_GREEN = 0, 1, 0, 0.5
HALF_BLUE  = 0, 0, 1, 0.5

# //////////////////////////////////////////////////////////////////////////////
