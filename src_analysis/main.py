import argparse
from _params import *
from calculations import *
from visualizations import *
from documentation import docs_main
from argparse import RawDescriptionHelpFormatter

# //////////////////////////////////////////////////////////////////////////////
def run_calculations():
    for run in RUNS:
        print(f"***** CURRENT RUN: {run} *****")
        PATH_GRO       = DIR_DA_TRAJECTORIES / f"{run}.gro"
        PATH_XTC       = DIR_DA_TRAJECTORIES / f"{run}.xtc"
        PATH_COORDS    = DIR_DA_COORDS / f"{run}.npy"
        PATH_COORDS_AS = DIR_DA_COORDS / f"{run}-AS0.npy"

        PATH_BSE_NAIVE  = DIR_DA_BSE  / f"{run}.naive.npy"
        PATH_BSE_ALTERN = DIR_DA_BSE  / f"{run}.altern.npy"
        PATH_CLUSTERING = DIR_DA_CLUSTERING / f"{run}.link.npy"
        PATH_PCA_SPACE  = DIR_DA_PCA  / f"{run}.space.npy"
        PATH_PCA_CUMVAR = DIR_DA_PCA  / f"{run}.cumvar.npy"
        PATH_PCA_PCOMPS = DIR_DA_PCA  / f"{run}.pcomponents.npy"
        PATH_RAMA       = DIR_DA_RAMA / f"{run}.npy"
        PATH_RGYR       = DIR_DA_RGYR / f"{run}.npy"
        PATH_RMSD       = DIR_DA_RMSD / f"{run}.npy"
        PATH_RMSF       = DIR_DA_RMSF / f"{run}.npy"


        traj = mda.Universe(str(PATH_GRO), str(PATH_XTC))
        coords = np.load(PATH_COORDS) if PATH_COORDS.exists() else extract_ca_coords(traj, PATH_COORDS)
        if not PATH_COORDS_AS.exists(): extract_AS_coords(traj, PATH_COORDS_AS, single_subunit = "t1" in run)

        if not PATH_RAMA.exists(): calc_rama(traj, PATH_RAMA)
        if not PATH_RGYR.exists(): calc_rgyr(traj, PATH_RGYR)
        if not PATH_RMSD.exists(): calc_rmsd(coords, PATH_RMSD)
        if not PATH_RMSF.exists(): calc_rmsf(traj, PATH_RMSF)

        if not PATH_PCA_SPACE.exists(): calc_pca(traj, PATH_PCA_SPACE, PATH_PCA_CUMVAR, PATH_PCA_PCOMPS)
        if not PATH_CLUSTERING.exists(): calc_link(PATH_RMSD, PATH_CLUSTERING, CLUSTERING_LINK_METHOD)

        if not PATH_BSE_NAIVE.exists(): calc_bse_naive(PATH_RMSD, PATH_BSE_NAIVE)
        if not PATH_BSE_ALTERN.exists(): calc_bse_alternate(PATH_RMSD, PATH_BSE_ALTERN)


# ------------------------------------------------------------------------------
def vis_bse():
    plotter = Plotter_BSE()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        plotter.vis_BSE(
            bse_naive_0 = np.load(DIR_DA_BSE / f"{run0}.naive.npy"),
            bse_altern_0 = np.load(DIR_DA_BSE / f"{run0}.altern.npy"),
            bse_naive_1 = np.load(DIR_DA_BSE / f"{run1}.naive.npy"),
            bse_altern_1 = np.load(DIR_DA_BSE / f"{run1}.altern.npy"),
            title = run_preffix,
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_cluster(plotter_objs):
    for run_preffix in RUN_PREFFIXES:
        plotter_objs.append(Plotter_RMSD1D_Clustering())
        run0 = run_preffix + "_rep0"

        plotter_objs[-1].vis_rmsd1d_clustering(
            rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy"),
            Z = np.load(DIR_DA_CLUSTERING / f"{run0}.link.npy"),
            color_map = COLOR_MAP_CLUSTERS,
            title = run0
        )


# ------------------------------------------------------------------------------
def vis_cmap0():
    plotter = Plotter_CMAP()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        coords0 = np.load(DIR_DA_COORDS / f"{run0}.npy")
        coords1 = np.load(DIR_DA_COORDS / f"{run1}.npy")

        plotter.base_cmap(coords0, frame = 0, title = f"{run0}-f0")
        plotter.base_cmap(coords1, frame = 0, title = f"{run1}-f0")


# ------------------------------------------------------------------------------
def vis_cmap1(plotter_objs):
    #### Example of specific analysis for these two trajectories
    coords_mt2 = np.load(DIR_DA_COORDS / f"mt2_rep0.npy")
    coords_wt2 = np.load(DIR_DA_COORDS / f"wt2_rep0.npy")

    # -------------------------------------------------------------------------- SPECIFIC ANALYSIS PART 1
    ### Example: Observe how CMAP values change in time. Not very insightful however
    plotter_objs.append(Plotter_CMAP_Dynamic())
    plotter_objs[-1].viscmap_interactive(coords_mt2, "mt2_rep0")

    ### Idea: Compare change of CMAP values from one frame relative to another reference frame
    plotter_objs.append(Plotter_DCMD_1Traj_2Frame())
    plotter_objs[-1].viscmap_interactive(coords_mt2, "mt2_rep0")

    ### Same as before but limiting the frames to a time section of interest
    plotter_objs.append(Plotter_DCMD_1Traj_2Frame())
    plotter_objs[-1].viscmap_interactive(coords_mt2, "mt2_rep0", min_frame = 4000, max_frame = 6000)

    ### Instead of comparing two frames from the same trajectory,
    ### it's also possible to compare the same frame in two different trajetories
    plotter_objs.append(Plotter_DCMD_2Traj_1Frame())
    plotter_objs[-1].viscmap_interactive(coords_mt2, coords_wt2, "mt2_rep0 vs wt2_rep0", min_frame = 4000, max_frame = 6000)

    # -------------------------------------------------------------------------- SPECIFIC ANALYSIS PART 2
    ### With DCMD, frames of interest can be isolated and further inspected
    frame_before = 5472
    frame_during = 5510
    frame_after = 5692

    ### Compare change of CMAP values between two frames (fixed)
    plotter = Plotter_CMAP()
    plotter.diff_cmap(coords_mt2, frame_before, frame_during, f"mt2_rep0 (f{frame_before} vs f{frame_during})")
    plotter.diff_cmap(coords_mt2, frame_during, frame_after, f"mt2_rep0 (f{frame_during} vs f{frame_after})")
    plotter.diff_cmap(coords_mt2, frame_before, frame_after, f"mt2_rep0 (f{frame_before} vs f{frame_after})")


# ------------------------------------------------------------------------------
def vis_cmap2():
    for rep in (0, 1):
        Plotter_CMAP_AS(
            coords_AS = {
                "mt1_as0" : np.load(DIR_DA_COORDS / f"mt1_rep{rep}-AS0.npy"),
                "mt2_as0" : np.load(DIR_DA_COORDS / f"mt2_rep{rep}-AS0.npy"),
                "mt2_as1" : np.load(DIR_DA_COORDS / f"mt2_rep{rep}-AS1.npy"),
                "wt1_as0" : np.load(DIR_DA_COORDS / f"wt1_rep{rep}-AS0.npy"),
                "wt2_as0" : np.load(DIR_DA_COORDS / f"wt2_rep{rep}-AS0.npy"),
                "wt2_as1" : np.load(DIR_DA_COORDS / f"wt2_rep{rep}-AS1.npy")
            },
            title = f"Average distance between interacting atoms of the Active Site (rep {rep})"
        )


# ------------------------------------------------------------------------------
def vis_pca(plotter_objs):
    for run_preffix in RUN_PREFFIXES:
        plotter_objs.append(Plotter_PCA())
        run0 = run_preffix + "_rep0"

        ### Alternatively, can use the 'vis_1pca_plotly' method
        ### With it, interaction with 3D scatter plot is faster
        ### However, it opens several windows, so it's advised to use it for only 1 run

        plotter_objs[-1].vis_1pca(
            cumvar = np.load(DIR_DA_PCA / f"{run0}.cumvar.npy"),
            space = np.load(DIR_DA_PCA / f"{run0}.space.npy"),
            pcomps = np.load(DIR_DA_PCA / f"{run0}.pcomponents.npy"),
            title = run_preffix
        )


# ------------------------------------------------------------------------------
def vis_pyinteraph():
    plotter = Plotter_Pyinteraph()
    # plotter.vis_pyinteraph(
    #     path_reduced = DIR_DA_PYINTERAPH / "hydrogen-bonds_all.csv",
    #     title = "Hydrogen Bonds"
    # )
    plotter.vis_pyinteraph(
        path_reduced = DIR_DA_PYINTERAPH / "salt-bridges_all.csv",
        title = "Salt Bridges"
    )
    plotter.vis_pyinteraph(
        path_reduced = DIR_DA_PYINTERAPH / "hydrophobic-clusters_all.csv",
        title = "Hydrophobic Clusters"
    )


# ------------------------------------------------------------------------------
def vis_rama(plotter_objs):
    for run_preffix in RUN_PREFFIXES:
        plotter_objs.append(Plotter_RAMA())
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        plotter_objs[-1].vis_2rama(
            rama0 = np.load(DIR_DA_RAMA / f"{run0}.npy"),
            rama1 = np.load(DIR_DA_RAMA / f"{run1}.npy"),
            title = run_preffix,
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_rgyr():
    plotter = Plotter_RGYR()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"
        plotter.vis_2rgyr(
            rgyr0 = np.load(DIR_DA_RGYR / f"{run0}.npy"),
            rgyr1 = np.load(DIR_DA_RGYR / f"{run1}.npy"),
            title = run_preffix,
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_rmsd0():
    """RMSD 2D (memory intensive)"""

    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy")
        rmsd_mat1 = np.load(DIR_DA_RMSD / f"{run1}.npy")

        plotter_rmsd2d = Plotter_RMSD2D()
        plotter_rmsd2d.vis_rmsd2d(rmsd_mat0, title = run0)
        plotter_rmsd2d.vis_rmsd2d(rmsd_mat1, title = run1)


# ------------------------------------------------------------------------------
def vis_rmsd1(plotter_objs):
    """RMSD 1D compare (intra comparisons)"""

    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        rmsd_mat0 = np.load(DIR_DA_RMSD / f"{run0}.npy")
        rmsd_mat1 = np.load(DIR_DA_RMSD / f"{run1}.npy")

        plotter_objs.append(Plotter_RMSD1D_Compare())
        plotter_objs[-1].vis_rmsd1d(
            rmsd_mat0, rmsd_mat1, run_preffix, run0, run1
        )


# ------------------------------------------------------------------------------
def vis_rmsf0():
    """RMSF 1D compare (intra comparisons)"""

    plotter = Plotter_RMSF()
    for run_preffix in RUN_PREFFIXES:
        run0 = run_preffix + "_rep0"
        run1 = run_preffix + "_rep1"

        plotter.vis_2rmsf(
            rmsf0 = np.load(DIR_DA_RMSF / f"{run0}.npy"),
            rmsf1 = np.load(DIR_DA_RMSF / f"{run1}.npy"),
            title = run0[:3] if (run0[:3] == run1[:3]) else f"{run0} vs {run1}",
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_rmsf1():
    """RMSF 1D compare (inter comparisons)"""

    comparisons = [
        ("mt1_rep0", "wt1_rep0"),
        ("mt2_rep0", "wt2_rep0"),
        ("mt1_rep0", "mt2_rep0"),
        ("wt1_rep0", "wt2_rep0"),
    ]

    plotter = Plotter_RMSF()
    for run0, run1 in comparisons:
        plotter.vis_2rmsf(
            rmsf0 = np.load(DIR_DA_RMSF / f"{run0}.npy"),
            rmsf1 = np.load(DIR_DA_RMSF / f"{run1}.npy"),
            title = run0[:3] if (run0[:3] == run1[:3]) else f"{run0} vs {run1}",
            label0 = run0,
            label1 = run1
        )


# ------------------------------------------------------------------------------
def vis_sasa():
    plotter = Plotter_SASA()
    plotter.vis_sasa(
        path_xvg = DIR_DA_SASA / "area.xvg",
        title = "Solvent Accessible Surface",
        xlabel = "Time ($ps$)",
        ylabel = r"Area ($nm^2$)"
    )
    plotter.vis_sasa_2vals(
        path_xvg = DIR_DA_SASA / "oa.xvg",
        title = "Area per atom over the trajectory",
        xlabel = "Atom",
        ylabel0 = "Average Area ($nm^2$)",
        ylabel1 = "Std Dev Area ($nm^2$)"
    )
    plotter.vis_sasa_2vals(
        path_xvg = DIR_DA_SASA / "or.xvg",
        title = "Area per residue over the trajectory",
        xlabel = "Residue",
        ylabel0 = "Average Area ($nm^2$)",
        ylabel1 = "Std Dev Area ($nm^2$)"
    )
    plotter.vis_sasa(
        path_xvg = DIR_DA_SASA / "sfe.xvg",
        title = "Free Energy of Solvation",
        xlabel = "Time ($ps$)",
        ylabel = "D Gsolv"
    )
    plotter.vis_sasa_2vals(
        path_xvg = DIR_DA_SASA / "volume.xvg",
        title = "Volume and Density",
        xlabel = "Time ($ps$)",
        ylabel0 = "Volume ($nm^3$)",
        ylabel1 = "Density ($g/l$)"
    )


# //////////////////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = docs_main["main"], formatter_class = RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--calc", action = "store_true", help = docs_main["calc"])
    parser.add_argument("-b", "--bse", action = "store_true", help = docs_main["bse"])
    parser.add_argument("-l", "--cluster", action = "store_true", help = docs_main["cluster"])
    parser.add_argument("-m0", "--cmap0", action = "store_true", help = docs_main["cmap0"]) # GENERAL ANALYSIS
    parser.add_argument("-m1", "--cmap1", action = "store_true", help = docs_main["cmap1"]) # SPECIFIC ANALYSIS
    parser.add_argument("-m2", "--cmap2", action = "store_true", help = docs_main["cmap2"]) # SPECIFIC ANALYSIS ACTIVE SITE
    parser.add_argument("-p", "--pca", action = "store_true", help = docs_main["pca"])
    parser.add_argument("-y", "--pyinteraph", action = "store_true", help = docs_main["pyinteraph"])
    parser.add_argument("-a", "--rama", action = "store_true", help = docs_main["rama"])
    parser.add_argument("-g", "--rgyr", action = "store_true", help = docs_main["rgyr"])
    parser.add_argument("-d0", "--rmsd0", action = "store_true", help = docs_main["rmsd0"]) # 2D
    parser.add_argument("-d1", "--rmsd1", action = "store_true", help = docs_main["rmsd1"]) # 1D
    parser.add_argument("-f0", "--rmsf0", action = "store_true", help = docs_main["rmsf0"]) # COMPARE SAME SYSTEM
    parser.add_argument("-f1", "--rmsf1", action = "store_true", help = docs_main["rmsf1"]) # COMPARE DIFFERENT SYSTEMS
    parser.add_argument("-s", "--sasa", action = "store_true", help = docs_main["sasa"])
    args = parser.parse_args()

    plot_flags = {**vars(args)}
    plot_flags.pop("calc")
    plot_requested = any(plot_flags.values())

    if not (args.calc | plot_requested):
        print(">>> Refer to the documentation (main.py -h) for description of usage. Closing...")
        exit()

    if args.calc: run_calculations()

    if plot_requested:
        if not any(DIR_DA_COORDS.iterdir()):
            print("xxx Visualization(s) requested, but it seems calculations haven't been perform yet (_coords folder empty). Aborting...")
            exit()
    else:
        print("/// No visualizations requested, closing program...")
        exit()

    ### must keep an active pointer of the plotter objects to avoid
    ### garbage collector to dispose of them, which would break interactivity
    ### i.e. functionality of widgets (sliders)
    plotter_objs = []
    if args.bse: vis_bse()
    if args.cluster: vis_cluster(plotter_objs)
    if args.cmap0: vis_cmap0() # GENERAL ANALYSIS
    if args.cmap1: vis_cmap1(plotter_objs) # SPECIFIC ANALYSIS
    if args.cmap2: vis_cmap2() # SPECIFIC ANALYSIS ACTIVE SITE
    if args.pca: vis_pca(plotter_objs)
    if args.pyinteraph: vis_pyinteraph()
    if args.rama: vis_rama(plotter_objs)
    if args.rgyr: vis_rgyr()
    if args.rmsd0: vis_rmsd0() # 2D
    if args.rmsd1: vis_rmsd1(plotter_objs) # 1D
    if args.rmsf0: vis_rmsf0() # COMPARE SAME SYSTEM
    if args.rmsf1: vis_rmsf1() # COMPARE DIFFERENT SYSTEMS
    if args.sasa: vis_sasa()

    plt.show()


# //////////////////////////////////////////////////////////////////////////////
