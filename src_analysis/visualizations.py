from _params import *
from calculations import calc_cluster, calc_cmap, calc_cmap_AS

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.widgets import Slider, Button

from MDAnalysis.analysis.data.filenames import Rama_ref

# //////////////////////////////////////////////////////////////////////////////
class Plotter:
    def stylize_ax(self, ax, title, xlabel, ylabel):
        ax.set_title(title, fontdict = dict(fontsize = 20))
        ax.set_xlabel(xlabel, fontdict = dict(fontsize = 16))
        ax.set_ylabel(ylabel, fontdict = dict(fontsize = 16))
        ax.tick_params(labelsize = 12)


    def plot_scatter(self, ax, x, y, colors, marker = '.'):
        line = ax.scatter(x, y, c = colors, marker = marker)
        return line


    def plot_heatmap(self, fig, ax, data, cmap):
        im = ax.imshow(data, cmap)
        colorbar = fig.colorbar(im)
        colorbar.ax.tick_params(labelsize = 12)
        return im, colorbar


    def plot_heatmap_cax(self, cax, ax, data, cmap):
        im = ax.imshow(data, cmap)
        colorbar = plt.colorbar(im, cax = cax)
        colorbar.ax.tick_params(labelsize = 12)
        return im, colorbar


    def update_line_data(self, line, x, data):
        line.set_offsets(
            np.append(x, data).reshape((2, x.size)).T
        )


# ////////////////////////////////////////////////////////////////////////////// BSE
class Plotter_BSE(Plotter):
    def vis_BSE(self, bse_naive_0, bse_altern_0, bse_naive_1, bse_altern_1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting BSE for '{title}'...")

        _, ax = plt.subplots()
        self.stylize_ax(ax, title, "Chunk Size", "Standard Deviation")
        ax.plot(bse_naive_0, color = BLUE, linestyle = ':', label = label0 + "_naive")
        ax.plot(bse_altern_0, color = BLUE, linestyle = '-', label = label0 + "_altern")
        ax.plot(bse_naive_1, color = RED, linestyle = ':', label = label1 + "_naive")
        ax.plot(bse_altern_1, color = RED, linestyle = '-', label = label1 + "_altern")
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// CMAP
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_CMAP(Plotter):
    do_abs_on_diff = True

    def __init__(self):
        self.do_abs = abs if self.do_abs_on_diff else (lambda x:x)


    def base_cmap(self, coords, frame, title = ''):
        print(f">>> Plotting base CMAP for '{title}'...")
        self.plot_cmap(title, calc_cmap(coords[frame]), "Reds_r")


    def diff_cmap(self, coords, frame0, frame1, title = '', color_method = "Blues"):
        print(f">>> Plotting diff CMAP for '{title}' (frame {frame1} - frame {frame0})...")
        cmap0 = calc_cmap(coords[frame0])
        cmap1 = calc_cmap(coords[frame1])
        self.plot_cmap(title, self.do_abs(cmap1 - cmap0), color_method)


    def plot_cmap(self, title, mat, color_method):
        fig, ax = plt.subplots()
        self.stylize_ax(ax, title, "CA Atom", "CA Atom")
        self.plot_heatmap(fig, ax, mat, color_method)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_CMAP_Dynamic(Plotter_CMAP):
    # -------------------------------------------------------------------------- OVERRIDES
    def base_cmap(self, frame):
        self.plot_cmap(calc_cmap(self.coords0[frame]), "Reds_r")


    def diff_cmap(self, frame0, frame1, color_method = "Blues"):
        cmap0 = calc_cmap(self.coords0[frame0])
        cmap1 = calc_cmap(self.coords0[frame1])
        self.plot_cmap(self.do_abs(cmap1 - cmap0), color_method)


    def plot_cmap(self, cmap, color_method):
        self.im, self.colorbar = self.plot_heatmap_cax(self.ax_dict['x'], self.ax_dict['a'], cmap, color_method)


    # --------------------------------------------------------------------------
    def viscmap_interactive(self, coords0, title = '', min_frame = 0, max_frame = 6000):
        print(f">>> Plotting base CMAP for '{title}'...")

        ##### DATA
        self.min_frame = min_frame
        self.max_frame = max_frame
        self.coords0 = coords0

        ##### FIGURES CREATION
        self.init_axes()
        self.stylize_ax(self.ax_dict['a'], title, "CA Atom", "CA Atom")

        ##### PLOTTING
        self.init_plots()

        ##### WIDGETS
        self.init_widgets()


    # --------------------------------------------------------------------------
    def init_axes(self):
        self.fig = plt.figure(layout = "constrained")
        self.ax_dict = self.fig.subplot_mosaic('\n'.join([
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "bbbbbbbbbb"
        ]))

    def init_plots(self):
        self.base_cmap(self.min_frame)

    def init_widgets(self):
        self.init_slider_frame0()

    # --------------------------------------------------------------------------
    def init_slider_frame0(self):
        self.slid_frame0 = Slider(
            ax = self.ax_dict['b'],
            label = "frame", color = "orange" ,
            valstep = 1, valinit = self.min_frame,
            valmin = self.min_frame, valmax = self.max_frame,
            valfmt = "%04i"
        )
        self.slid_frame0.on_changed(self.update_cmap)

    def update_cmap(self, val):
        frame = self.slid_frame0.val
        cmap = calc_cmap(self.coords0[frame])

        self.im.set_data(cmap)
        self.im.set_clim(vmin = 0, vmax = np.max(cmap))
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_DCMD_1Traj_2Frame(Plotter_CMAP_Dynamic):
    def init_axes(self):
        self.fig = plt.figure(layout = "constrained")
        self.ax_dict = self.fig.subplot_mosaic('\n'.join([
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "aaaaaaaaax",
            "bbbbbbbbbb",
            "cccccccccc"
        ]))

    def init_plots(self):
        self.diff_cmap(self.min_frame, self.min_frame + 1, "seismic")

    def init_widgets(self):
        self.init_slider_frame0()
        self.init_slider_frame1()

    # --------------------------------------------------------------------------
    def init_slider_frame1(self):
        self.slid_frame1 = Slider(
            ax = self.ax_dict['c'],
            label = "frame", color = "blue",
            valstep = 1, valinit = self.min_frame + 1,
            valmin = self.min_frame, valmax = self.max_frame,
            valfmt = "%04i"
        )
        self.slid_frame1.on_changed(self.update_cmap)

    def update_cmap(self, val):
        frame0 = self.slid_frame0.val
        frame1 = self.slid_frame1.val
        diff = self.do_abs(calc_cmap(self.coords0[frame1]) - calc_cmap(self.coords0[frame0]))

        self.im.set_data(diff)
        self.im.set_clim(vmin = np.min(diff), vmax = np.max(diff))
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_DCMD_2Traj_1Frame(Plotter_CMAP_Dynamic):
    def viscmap_interactive(self, coords0, coords1, title = '', min_frame = 0, max_frame = 6000):
        self.coords1 = coords1
        super().viscmap_interactive(coords0, title, min_frame, max_frame)


    def init_plots(self):
        cmap0 = calc_cmap(self.coords0[self.min_frame])
        cmap1 = calc_cmap(self.coords1[self.min_frame])
        self.plot_cmap(self.do_abs(cmap1 - cmap0), "seismic")

    # --------------------------------------------------------------------------
    def update_cmap(self, frame):
        diff = self.do_abs(calc_cmap(self.coords1[frame]) - calc_cmap(self.coords0[frame]))

        self.im.set_data(diff)
        self.im.set_clim(0, vmax = np.max(diff))
        self.colorbar.update_normal(self.im)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_CMAP_AS(Plotter):
    def __init__(self, coords_AS, title = ''):
        sns.set_palette("bwr")
        fig, ax = plt.subplots()

        interactions = ["AB", "CD", "EF", "FG", "EG", "GH", "IJ"]
        df = pd.DataFrame(columns = ["Trajectory", "Interaction", "Distance"])

        print(f">>> Calculating CMAPs for '{title}'...")
        for system,coords in coords_AS.items():
            atoms = {atom : coords[:,i,:] for i,atom in enumerate("ABCDEFGHIJ")}
            for atom0, atom1 in interactions:
                new_df = pd.DataFrame()
                new_df["Distance"] = np.diag(calc_cmap_AS(atoms[atom0], atoms[atom1]))
                new_df["Interaction"] = atom0 + atom1
                new_df["Trajectory"] = system
                df = df.append(new_df)

        print(f"...>>> Plotting...")
        # sns.boxplot(data = df, x = "Interaction", y = "Distance", hue = "Trajectory", ax = ax)
        sns.pointplot(data = df, x = "Interaction", y = "Distance", hue = "Trajectory", ax = ax)
        plt.title(title)


# ////////////////////////////////////////////////////////////////////////////// PCA
class Plotter_PCA(Plotter):
    def vis_1pca(self, cumvar, space, pcomps, title = ''):
        print(f">>> Plotting PCA for '{title}'...")

        ##### DATA
        self.x = np.arange(pcomps.shape[0])
        self.pcomps = np.copy(pcomps)
        n_pcs = np.where(cumvar > 0.95)[0][0]

        print("...>>> Components:", len(cumvar))
        print("...>>> with cumvar > 0.95:", n_pcs)

        ##### FIGURES CREATION
        _, ax_cumvar = plt.subplots()
        self.stylize_ax(ax_cumvar, title, "Components", "Cumulative Variance")
        ax_cumvar.set_xlim([0, 100])

        fig = plt.figure()
        fleft, self.fright = fig.subfigures(ncols = 2)

        ax_3d = fleft.add_subplot(projection = "3d")
        ax_3d.set_title(title, fontdict = dict(fontsize = 20))
        self.ax_dict = self.fright.subplot_mosaic("a;a;a;b")
        self.stylize_ax(self.ax_dict['a'], title, "Component", "Variance")

        ##### WIDGETS
        self.slid_pc = Slider(
            ax = self.ax_dict['b'],
            label = "pc", color = "orange",
            valstep = 1, valinit = 0,
            valmin = 0, valmax = pcomps.shape[0] - 1
        )
        self.slid_pc.on_changed(self.update_plot)

        ##### PLOTTING
        ax_cumvar.plot(cumvar)
        self.line_2d = self.ax_dict['a'].scatter(self.x, pcomps[0], marker = '.')

        sc = ax_3d.scatter(space[:,0], space[:,1], space[:,2], c = np.arange(space.shape[0]), s = 20, alpha = 0.5)
        ax_3d.legend(*sc.legend_elements(), loc = "lower left")


    def update_plot(self, idx):
        self.update_line_data(self.line_2d, self.x, self.pcomps[idx])


    def vis_1pca_plotly(self, cumvar, space, pcomps, title = ''):
        import pandas as pd
        import plotly.express as px

        ##### CUMVAR PLOT
        fg = px.line(
            x = np.arange(cumvar.shape[0]),
            y = cumvar,
            labels = {"x" : "components", "y" : "cumulated variance"},
            range_x = [0, 100]
        )
        fg.show()

        ##### 3D SCATTER
        pca_data = pd.DataFrame(space, columns = ["first_comp", "second_comp", "third_comp"])
        fig = px.scatter_3d(
            pca_data, x = "first_comp", y = "second_comp", z = "third_comp",
            color = pca_data.index.values, width = 900, height = 800
        )
        fig.show()

        ##### 2D PLOTS
        fg0 = px.line(x = np.arange(pcomps.shape[0]), y = pcomps[0])
        fg0.show()

        fg1 = px.line(x = np.arange(pcomps.shape[0]), y = pcomps[1])
        fg1.show()

        fg2 = px.line(x = np.arange(pcomps.shape[0]), y = pcomps[2])
        fg2.show()


# ////////////////////////////////////////////////////////////////////////////// PYINTERAPH
class Plotter_Pyinteraph(Plotter):
    def vis_pyinteraph(self, path_reduced, title = ''):
        print(f">>> Plotting Pyinteraph for '{title}'...")

        ##### DATA
        header = ["res0_chain", "res0_num", "res0_name", "res0_inter_group", "res1_chain", "res1_num", "res1_name", "res1_inter_group", "occurrence"]
        df = pd.read_csv(path_reduced, header = None, usecols = range(len(header)), names = header)
        matrix = np.zeros((N_RESIDS_SUBUNIT, N_RESIDS_SUBUNIT))

        for _,row in df.iterrows():
            id0 = row.res0_num - 1
            id1 = row.res1_num - 1
            matrix[id0, id1] = row.occurrence
            matrix[id1, id0] = row.occurrence

        ##### FIGURES CREATION
        fig, ax = plt.subplots()
        self.stylize_ax(ax, title, "Residue", "Residue")

        ##### PLOTTING
        self.plot_heatmap(fig, ax, matrix, "hot")


# ////////////////////////////////////////////////////////////////////////////// RAMA
class Plotter_RAMA(Plotter):
    def vis_2rama(self, rama0, rama1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting RAMA for '{title}'...")

        ##### DATA
        self.rama_mat0 = rama0
        self.rama_mat1 = rama1

        ##### FIGURES CREATION
        self.init_axes()
        self.stylize_ax(self.ax_dict['a'], title, r"$\phi$", r"$\psi$")

        ##### PLOTTING
        self.init_plots(label0, label1)

        ##### WIDGETS
        self.init_widgets()


    def init_axes(self):
        fig = plt.figure(layout = "constrained")
        self.ax_dict = fig.subplot_mosaic("a;a;a;a;a;a;b")

        #### Rama ax formatting, adapted from "MDAnalysis/analysis/dihedrals.Ramachandran.plot"
        self.ax_dict['a'].set_facecolor((0,0,0,1))
        self.ax_dict['a'].axis([-180, 180, -180, 180])
        self.ax_dict['a'].axhline(0, color = 'w', lw = 1)
        self.ax_dict['a'].axvline(0, color = 'w', lw = 1)
        self.ax_dict['a'].set(xticks = range(-180, 181, 60), yticks = range(-180, 181, 60))

        degree_formatter = plt.matplotlib.ticker.StrMethodFormatter("{x:g}°")
        self.ax_dict['a'].xaxis.set_major_formatter(degree_formatter)
        self.ax_dict['a'].yaxis.set_major_formatter(degree_formatter)


    def init_plots(self, label0, label1):
        #### Reference Rama plot, adapted from "MDAnalysis/analysis/dihedrals.Ramachandran.plot"
        X, Y = np.meshgrid(
            np.arange(-180, 180, 4),
            np.arange(-180, 180, 4)
        )
        levels = [1, 17, 15000]
        colors = ['#888888', '#AAAAAA']
        self.ax_dict['a'].contourf(X, Y, np.load(Rama_ref), levels = levels, colors = colors)


        #### Observed Rama values
        rama_arr0 = self.rama_mat0[0].reshape(self.rama_mat0.shape[1], 2)
        rama_arr1 = self.rama_mat1[0].reshape(self.rama_mat1.shape[1], 2)

        self.line_rama0 = self.ax_dict['a'].scatter(
            rama_arr0[:,0], rama_arr0[:,1], color = HALF_RED,  marker = 'x', s = 20, label = label0
        )
        self.line_rama1 = self.ax_dict['a'].scatter(
            rama_arr1[:,0], rama_arr1[:,1], color = HALF_BLUE, marker = '+', s = 20, label = label1
        )
        self.ax_dict['a'].legend(loc = "upper right")


    def init_widgets(self):
        self.slid_ref_frame = Slider(
            ax = self.ax_dict['b'],
            label = "ref_frame", color = "orange",
            valstep = 1, valinit = 0,
            valmin = 0, valmax = self.rama_mat0.shape[0] - 1,
            valfmt = "%04i"
        )
        self.slid_ref_frame.on_changed(self.update_plot)


    def update_plot(self, ref_frame):
        rama_arr0 = self.rama_mat0[ref_frame]
        rama_arr1 = self.rama_mat1[ref_frame]
        self.update_line_data(self.line_rama0, rama_arr0[:,0], rama_arr0[:,1])
        self.update_line_data(self.line_rama1, rama_arr1[:,0], rama_arr1[:,1])


# ////////////////////////////////////////////////////////////////////////////// RGYR
class Plotter_RGYR(Plotter):
    def vis_2rgyr(self, rgyr0, rgyr1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting RGYR for '{title}'...")

        frames = np.arange(rgyr0.size)
        _, ax = plt.subplots()
        self.stylize_ax(ax, title, "Frame", "Radius of Gyration")
        ax.scatter(frames, rgyr0, color = BLUE, marker = '+', s = 20, label = label0)
        ax.scatter(frames, rgyr1, color = HALF_RED, marker = 'x', s = 20, label = label1)
        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// RMSD
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_RMSD2D(Plotter):
    def vis_rmsd2d(self, rmsd_mat0, title = ''):
        print(f">>> Plotting RMSD 2D for '{title}'...")

        fig, ax = plt.subplots()
        self.stylize_ax(ax, title, "Frame", "Frame")
        self.plot_heatmap(fig, ax, rmsd_mat0, cmap = "Reds")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_RMSD1D(Plotter):
    init_ref_frame = 0

    def vis_rmsd1d(self, rmsd_mat0, title = ''):
        print(f">>> Plotting RMSD 1D for '{title}'...")

        ##### DATA
        self.x = np.arange(rmsd_mat0.shape[0])
        self.rmsd_mat0 = np.copy(rmsd_mat0)
        self.colors = np.zeros((rmsd_mat0.shape[0], 4))
        self.update_values(self.init_ref_frame)

        ##### FIGURES CREATION
        self.init_axes()
        self.stylize_ax(self.ax_dict['a'], title, "Frame", "RMSD")

        ##### PLOTTING
        self.init_plots()

        ##### WIDGETS
        self.init_widgets()

    # --------------------------------------------------------------------------
    def update_values(self, ref_frame):
        self.rmsd_arr0 = self.rmsd_mat0[ref_frame]
        self.colors[:] = BLUE
        self.colors[ref_frame] = RED

    def init_axes(self):
        fig = plt.figure(layout = "constrained")
        self.ax_dict = fig.subplot_mosaic("a;a;a;a;a;a;b")

    def init_plots(self):
        self.line_rmsd0 = self.ax_dict['a'].scatter(self.x, self.rmsd_arr0, color = self.colors, marker = '.')

    def init_widgets(self):
        self.init_slider_ref_frame()

    def init_slider_ref_frame(self):
        self.slid_ref_frame = Slider(
            ax = self.ax_dict['b'],
            label = "ref_frame", color = "orange",
            valstep = 1, valinit = self.init_ref_frame,
            valmin = 0, valmax = self.rmsd_mat0.shape[0] - 1,
            valfmt = "%04i"
        )
        self.slid_ref_frame.on_changed(self.update_plot)

    # --------------------------------------------------------------------------
    def update_plot(self, ref_frame):
        self.update_values(ref_frame)
        self.update_line_rmsd()

    def update_line_rmsd(self):
        self.update_line_data(self.line_rmsd0, self.x, self.rmsd_arr0)
        self.line_rmsd0.set_color(self.colors)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class Plotter_RMSD1D_Compare(Plotter_RMSD1D):
    def vis_rmsd1d(self, rmsd_mat0, rmsd_mat1, title = '', label0 = '', label1 = ''):
        self.rmsd_mat1 = rmsd_mat1
        self.label0 = label0
        self.label1 = label1
        super().vis_rmsd1d(rmsd_mat0, title)
        print(f"...>>> Comparison mode.")

    # --------------------------------------------------------------------------
    def update_values(self, ref_frame):
        super().update_values(ref_frame)
        self.rmsd_arr1 = self.rmsd_mat1[ref_frame]

    def init_plots(self):
        self.line_rmsd0 = self.ax_dict['a'].scatter(self.x, self.rmsd_arr0, color = BLUE, marker = '+', s = 12, label = self.label0)
        self.line_rmsd1 = self.ax_dict['a'].scatter(self.x, self.rmsd_arr1, color = HALF_RED, marker = 'x', s = 12, label = self.label1)
        self.ax_dict['a'].legend()

    def update_line_rmsd(self):
        self.update_line_data(self.line_rmsd0, self.x, self.rmsd_arr0)
        self.update_line_data(self.line_rmsd1, self.x, self.rmsd_arr1)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CLUSTERING
class Plotter_RMSD1D_Clustering(Plotter_RMSD1D):
    init_clustering_t = 900
    label_criterion = "distance"

    def vis_rmsd1d_clustering(self, rmsd_mat0, Z, color_map, title = ''):
        self.Z = Z
        self.color_map = color_map
        self.highlight_cluster = False
        self.update_values_cluster(self.init_clustering_t)

        super().vis_rmsd1d(rmsd_mat0, title)
        print(f"...>>> Clustering mode.")


    def update_values_cluster(self, t):
        self.cluster_mat = calc_cluster(self.Z, t, label_criterion = self.label_criterion)
        self.n_clusters = np.max(self.cluster_mat)
        self.cluster_colors = cm.__dict__[self.color_map](np.linspace(0, 1, self.n_clusters))

    # --------------------------------------------------------------------------
    def init_axes(self):
        self.fig = plt.figure(layout = "constrained")
        self.ax_dict = self.fig.subplot_mosaic("aaaa;aaaa;aaaa;aaaa;aaaa;aaaa;bbbd;cccd")

    def init_widgets(self):
        super().init_widgets()

        self.slid_clustering_t = Slider(
            ax = self.ax_dict['c'],
            label = "t value", color = "blue",
            valstep = 25, valinit = self.init_clustering_t,
            valmin = 25, valmax = 2000,
            valfmt = "%04i"
        )
        self.slid_clustering_t.on_changed(self.update_clustering_t)

        self.button_toggle_alpha = Button(
            ax = self.ax_dict['d'],
            label = f"{self.n_clusters} clusters"
        )
        self.button_toggle_alpha.on_clicked(self.toggle_alpha)

    # --------------------------------------------------------------------------
    def update_clustering_t(self, clustering_t):
        self.update_values_cluster(clustering_t)
        self.button_toggle_alpha.label.set_text(f"{self.n_clusters} clusters")
        self.update_plot(self.slid_ref_frame.val)

    def toggle_alpha(self, event):
        self.highlight_cluster = not self.highlight_cluster
        self.ax_dict['a'].set_facecolor((0, 0, 0) if self.highlight_cluster else (1, 1, 1))
        self.fig.canvas.draw_idle()
        self.update_plot(self.slid_ref_frame.val)

    def update_values(self, ref_frame):
        self.rmsd_arr0 = self.rmsd_mat0[ref_frame]
        current_cluster = self.cluster_mat[ref_frame]

        for i,color in enumerate(self.cluster_colors[:,:]):
            r,g,b,a = color
            if (i + 1 != current_cluster) and self.highlight_cluster:
                a = .1
            self.colors[self.cluster_mat == i + 1] = r,g,b,a

    # --------------------------------------------------------------------------


# ////////////////////////////////////////////////////////////////////////////// RMSF
class Plotter_RMSF(Plotter):
    def vis_2rmsf(self, rmsf0, rmsf1, title = '', label0 = '', label1 = ''):
        print(f">>> Plotting RMSF for '{title}'...")

        frames = np.arange(max(rmsf0.size, rmsf1.size))
        _, ax = plt.subplots()
        self.stylize_ax(ax, title, "CA atom", "RMSF")

        if rmsf0.size == rmsf1.size:
            ax.plot(frames, rmsf0, color = BLUE, linestyle = ':', linewidth = 2, label = label0)
            ax.plot(frames, rmsf1, color = HALF_RED, linestyle = '-', linewidth = 2, label = label1)
        else:
            rmsf_half, rmsf_full = sorted([rmsf0, rmsf1], key = lambda arr: arr.size)
            label_half, label_full = (label0, label1) if rmsf0 is rmsf_half else (label1, label0)

            ax.plot(frames, rmsf_full, color = BLUE, linestyle = ':', linewidth = 2, label = label_full)
            ax.plot(frames[:rmsf_half.size], rmsf_half, color = HALF_RED, linestyle = '-', linewidth = 2, label = label_half)
            ax.plot(frames[rmsf_half.size:], rmsf_half, color = HALF_GREEN, linestyle = '-', linewidth = 2, label = label_half)

        ax.legend()


# ////////////////////////////////////////////////////////////////////////////// SASA
class Plotter_SASA(Plotter):
    def vis_sasa(self, path_xvg, title = '', xlabel = '', ylabel = ''):
        print(f">>> Plotting SASA for '{title}'...")

        _, ax = plt.subplots()
        self.stylize_ax(ax, title, xlabel, ylabel)

        x,y = np.loadtxt(path_xvg, comments = ['#', '@'], unpack = True)
        ax.scatter(x, y, color = HALF_BLUE, marker = '+')


    def vis_sasa_2vals(self, path_xvg, title = '', xlabel = '', ylabel0 = '', ylabel1 = ''):
        print(f">>> Plotting SASA for '{title}'...")

        fig = plt.figure()
        ax_dict = fig.subplot_mosaic("ab")
        self.stylize_ax(ax_dict['a'], title, xlabel, ylabel0)
        self.stylize_ax(ax_dict['b'], title, xlabel, ylabel1)

        x,y0,y1 = np.loadtxt(path_xvg, comments = ['#', '@'], unpack = True)
        ax_dict['a'].scatter(x, y0, color = HALF_BLUE, marker = '+')
        ax_dict['b'].scatter(x, y1, color = HALF_RED,  marker = '+')


# //////////////////////////////////////////////////////////////////////////////
