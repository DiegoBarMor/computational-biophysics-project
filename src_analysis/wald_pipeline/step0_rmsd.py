from _params import *
from wald_pipeline._params import *
from wald_pipeline.info import Info
from visualizations import Plotter_RMSD1D

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# //////////////////////////////////////////////////////////////////////////////
class Plotter_RMSD1D_WALD(Plotter_RMSD1D):
    init_rmsd_treshold = 2

    def vis_rmsd1d(self, rmsd_mat0, title = ''):
        self.y = np.zeros(rmsd_mat0.shape[0])
        super().vis_rmsd1d(rmsd_mat0, title)
        print(f"...>>> WALD frame selection mode.")

    def init_axes(self):
        self.fig = plt.figure(layout = "constrained")
        self.ax_dict = self.fig.subplot_mosaic("aaaa;aaaa;aaaa;aaaa;aaaa;aaaa;bbbd;cccd")

    def init_plots(self):
        super().init_plots()
        self.line_treshold = self.ax_dict['a'].plot(self.x, self.y)

    def init_widgets(self):
        super().init_widgets()

        self.slid_rmsd_treshold = Slider(
            ax = self.ax_dict['c'],
            label = "rmsd_treshold", color = "blue",
            valstep = .05, valinit = self.init_rmsd_treshold,
            valmin = 0, valmax = np.max(self.rmsd_mat0),
            valfmt = "%02.2f"
        )
        self.slid_rmsd_treshold.on_changed(self.update_plot)

        self.button_save = Button(
            ax = self.ax_dict['d'],
            label = "Save WA frames"
        )
        self.button_save.on_clicked(self.save_selected_frames)


    def update_plot(self, val):
        self.update_values(self.slid_ref_frame.val, self.slid_rmsd_treshold.val)
        self.update_line_rmsd()
        self.update_line_treshold()

    def update_values(self, ref_frame, rmsd_treshold = None):
        if rmsd_treshold is None:
            rmsd_treshold = self.init_rmsd_treshold

        self.rmsd_arr0 = self.rmsd_mat0[ref_frame]

        self.colors[:] = BLUE
        self.colors[self.rmsd_arr0 <= rmsd_treshold] = GREEN
        self.colors[ref_frame] = RED

        self.y[:] = rmsd_treshold

    # --------------------------------------------------------------------------
    def update_line_treshold(self):
        self.line_treshold[0].set_data(self.x, self.y)

    # --------------------------------------------------------------------------
    def save_selected_frames(self, event):
        ref_frame = self.slid_ref_frame.val
        rmsd_treshold = self.slid_rmsd_treshold.val

        self.rmsd_arr0 = self.rmsd_mat0[ref_frame]
        frames = [int(f) for f in self.x[self.rmsd_arr0 <= rmsd_treshold]]

        info = Info(DIR_DA_WALD / f"{RUN_WALD}-info.json")

        info.update(
            rsmd_treshold = rmsd_treshold,
            ref_frame = ref_frame,
            ref_frame_index = frames.index(ref_frame),
            frames = frames,
        )

        self.button_save.label.set_text(f"Saved {len(frames)} frames")


# //////////////////////////////////////////////////////////////////////////////
