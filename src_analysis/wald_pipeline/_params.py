from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////
RUN_WALD = "mt2_rep0"
WALD_DIVISIONS = 20
WALD_DENSITY_TRESHOLD_LOWER = 0
WALD_DENSITY_TRESHOLD_UPPER = 3

###########################################################
DIR_DA = Path.cwd().parent / "data_analysis"
DIR_DA_TRAJECTORIES = DIR_DA / "_trajectories"
DIR_DA_WALD = DIR_DA / "wald"

WALD_NAME = RUN_WALD
PATH_WALD_INFO = DIR_DA_WALD / f"{WALD_NAME}-info.json"


# //////////////////////////////////////////////////////////////////////////////
