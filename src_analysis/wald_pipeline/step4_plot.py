from wald_pipeline._params import *

import numpy as np
import pandas as pd
import plotly.graph_objects as go

# //////////////////////////////////////////////////////////////////////////////
def load_vertices_faces(info):
    vertices = []; faces = []
    for i in range(info.n_repr):
        vertices_dir = DIR_DA_WALD / f"{WALD_NAME}-mesh_v{i}.npy"
        faces_dir =    DIR_DA_WALD / f"{WALD_NAME}-mesh_f{i}.csv"

        vertices.append(np.load(vertices_dir))
        faces.append(pd.read_csv(faces_dir))

        vertices[i][:,0] *= info.x_scale
        vertices[i][:,1] *= info.y_scale
        vertices[i][:,2] *= info.z_scale

        vertices[i][:,0] -= info.x_offset
        vertices[i][:,1] -= info.y_offset
        vertices[i][:,2] -= info.z_offset

    return vertices, faces

def plot_layer(fig, vertices, faces = None):
    kwargs = dict(
        x = vertices[:,0],
        y = vertices[:,1],
        z = vertices[:,2],
    )

    if faces is not None:
        kwargs['i'] = faces.i
        kwargs['j'] = faces.j
        kwargs['k'] = faces.k
        kwargs["facecolor"] = faces.color

    fig.add_trace(
        go.Mesh3d(**kwargs)
    )

def plot_density(fig, info):
    density = np.load(DIR_DA_WALD / f"{WALD_NAME}-wald.npy")

    X, Y, Z = np.mgrid[
        0 : info.xbox : complex(0, WALD_DIVISIONS),
        0 : info.ybox : complex(0, WALD_DIVISIONS),
        0 : info.zbox : complex(0, WALD_DIVISIONS),
    ]

    fig.add_trace(
        go.Volume(
            x = X.flatten(),
            y = Y.flatten(),
            z = Z.flatten(),
            value =  -density.flatten(),
            isomin = -WALD_DENSITY_TRESHOLD_UPPER,
            isomax = -WALD_DENSITY_TRESHOLD_LOWER,
            opacity = 0.05, # needs to be small to see through all surfaces
            opacityscale = "uniform",
            surface_count = 50, # needs to be a large number for good volume rendering
            colorscale = ["blue", "white", "magenta"],
        )
    )

def plot_protein(fig, info):
    for v,f in zip(*load_vertices_faces(info)):
        plot_layer(fig, v, f)


# //////////////////////////////////////////////////////////////////////////////
