import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
def gen_cut_xtc(gro_path, xtc_path, wald_gro_path, wald_xtc_path, info):
    print(">>> Cutting the frames of interest...")

    traj = mda.Universe(str(gro_path), str(xtc_path))
    protein = traj.select_atoms("protein")

    with mda.Writer(str(wald_gro_path), protein.n_atoms) as W:
        W.write(protein)

    with mda.Writer(str(wald_xtc_path), protein.n_atoms) as W:
        for ts in traj.trajectory[info.frames]:
            W.write(protein)


def gen_vmd(gro_path, xtc_path, vmd_path, info):
    print(">>> Generating VMD visualization state...")

    gro_path = str(gro_path).replace('\\', '/')
    xtc_path = str(xtc_path).replace('\\', '/')

    out = """color change rgb red 0 1 0;   # acidic (hydrophilic)
color change rgb blue 0 1 0;  # basic (hydrophilic)
color change rgb green 0 1 0; # polar (hydrophilic)
color change rgb white 1 1 0; # non-polar (hydrophobic)
color change rgb cyan 1 0 0;  # active site
\n"""

    out += f"mol new {gro_path} type gro first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
    out += f"mol addfile {xtc_path} type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"
    out += """
mol delrep 0 top
mol representation QuickSurf
mol color ResType
mol selection {protein and not resid 41 145 164}
mol material Opaque
mol addrep top
mol smoothrep top 0 4

mol representation QuickSurf
mol color Element
mol selection {resid 41 145 164}
mol material Opaque
mol addrep top
mol smoothrep top 0 4

pbc box\n"""
    out += f"puts \">>> Reference frame at index {info.ref_frame_index}\""

    with open(vmd_path, 'w') as file: file.write(out)


# //////////////////////////////////////////////////////////////////////////////
