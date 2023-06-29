from wald_pipeline._params import *

import re
import numpy as np
import pandas as pd

# //////////////////////////////////////////////////////////////////////////////
####################################### MATERIALS
def extract_colors(mtl_path):
    with open(mtl_path, 'r') as file: mtl = file.read()
    print("...>>> Extracting colors...")

    mat_names_raw = re.findall("\nnewmtl.*", mtl)
    mat_colors_raw = re.findall("\nKd.*", mtl)

    mat_names = map(lambda s: s[8:], mat_names_raw)
    mat_colors = map(lambda s: f"rgb({ s[4:].replace(' ', ',') })", mat_colors_raw)
    colors = dict(zip(mat_names, mat_colors))

    return mat_names, mat_colors, colors

####################################### VERTICES
def extract_vertices(repr, out_name = ''):
    vertices_raw = re.findall("\nv\s.*", repr)
    vertices = np.array(list(map(lambda s: s[3:].split(' '), vertices_raw)), dtype = float)

    if out_name:
        print(f"...>>> Saving vertices to '{out_name}'...")
        np.save(DIR_DA_WALD / out_name, vertices, allow_pickle = False)
        return vertices.shape[0]

    return vertices

####################################### FACES
def prepare_faces(repr):
    faces_start = re.search("\nusemtl.*", repr)
    if faces_start is None: return
    return repr[faces_start.start():]

def extract_faces(faces_str, n_vertices, colors, out_name):
    re_materials = re.finditer("usemtl.*", faces_str)
    materials = map(lambda m: (m.start(), m.group()), re_materials)


    matches = np.array(list(re.finditer("[\s\-]\d+[\/]", faces_str)))
    matches = np.reshape(matches, (matches.size // 3, 3))

    faces = pd.DataFrame(
        columns = ['i', 'j', 'k'],
        data = np.vectorize( lambda m: n_vertices + int( m.group()[:-1] ) )(matches)
    )

    raw_start = np.vectorize( lambda m: m.start() )(matches[:,0])
    color = np.full(raw_start.size, "rgb(_.--,_.--,_.--)")

    for id,m in materials:
        color[raw_start > id] = colors[m[7:]]
    faces["color"] = color

    print(f"...>>> Saving faces to '{out_name}'...")
    faces.to_csv(DIR_DA_WALD / out_name, index = False)


# //////////////////////////////////////////////////////////////////////////////
