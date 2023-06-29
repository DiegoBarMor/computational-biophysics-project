import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances
from MDAnalysis.analysis.rms import rmsd, RMSF
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.analysis.pca import PCA, cosine_content
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

# //////////////////////////////////////////////////////////////////////////////
# ------------------------------------------------------------------------------ COORDINATES EXTRACTION
def extract_ca_coords(traj, path_coords):
    print(f">>> Preparing '_coords/{path_coords.name}'...")
    coords = np.array([traj.select_atoms("protein and name CA").positions for _ in traj.trajectory[:]])
    np.save(path_coords, coords, allow_pickle = False)
    return coords

def extract_AS_coords(traj, path_coords_as, single_subunit = False):
    """Extract Active Site coordinates specifically for 8DFN and 7SI9"""
    print(f">>> Preparing '_coords/{path_coords_as.name}' (Active Site)...")

    selections = [
        "resid 145 and name SG",
        "resid 41 and name NE2",
        "resid 41 and name CE1",
        "resid 164 and name O",
        "resid 41 and name ND1",
        "resid 164 and name ND?",
        "resid 187 and name OD2",
        "resid 40 and name NE",
        "resid 187 and name OD1",
        "resid 40 and name NH2",
    ]

    # create an atom group
    atoms = sum(traj.select_atoms(sel) for sel in selections)

    coords = np.zeros((len(traj.trajectory), len(atoms), 3))
    for i,_ in enumerate(traj.trajectory[:]):
        for j,atom in enumerate(atoms):
            coords[i,j] = atom.position

    if single_subunit:
        np.save(path_coords_as, coords, allow_pickle = False)
    else:
        np.save(path_coords_as, coords[:, 0::2, :], allow_pickle = False)
        np.save(path_coords_as.parent / f"{path_coords_as.stem}-AS1.npy", coords[:, 1::2, :], allow_pickle = False)


# ------------------------------------------------------------------------------ BSE
def calc_bse_naive(path_rmsd, path_bse):
    """Naive version of BSE: the mean value is not stable since we are leaving out the last block.
    Adapted from the script provided by professor Tubiana."""

    print(f">>> Preparing 'bse/{path_bse.name}'...")

    v = np.load(path_rmsd)[0]
    n = len(v)
    l = []
    for i in range(1,n//3):
        d = n//i
        w = v[:d*i].reshape(d,i)
        m = w.mean(axis=1)
        l.append(m.var()/(m.shape[0]-1))
    l = np.asarray(l)
    l = np.sqrt(l)

    np.save(path_bse, l, allow_pickle = False)
    print("...>>> Done.")

def calc_bse_alternate(path_rmsd, path_bse):
    """Returns a vector s(m) with the std-dev computed according to the Block Analysis method. The average is stable, needs more testing.
    Adapted from the script provided by professor Tubiana."""

    print(f">>> Preparing 'bse/{path_bse.name}'...")

    v = np.load(path_rmsd)[0]
    n = len(v)
    a = []
    s = []
    for i in range(1,n//3):
        d = n//i
        r = n-d*i
        w = v[:d*i].reshape(d,i)
        m = w.mean(axis=1)
        avg = m.mean()
        weights = np.full(d,i)
        if d*i != n:
            m_d = v[d*i:].mean()
            avg = (d*i*avg+r*m_d)/(d*i+r)
            weights=np.append(weights,r)
            m=np.append(m,m_d)
        m2 = (m - avg)**2
        weights = weights/weights.sum()
        m2 *= weights
        var = m2.sum()/(m2.shape[0]-1)
        s.append(var)
        a.append(avg)
    s = np.asarray(s)
    s = np.sqrt(s)
    a = np.asarray(a)

    # a is the average, check that it is constant.
    a -= np.mean(a)
    print(f"...>>> BSE_Alternate: Average fluctuates by {np.max(np.abs(a)):.2e}.")

    np.save(path_bse, s, allow_pickle = False)
    print("...>>> Done.")


# ------------------------------------------------------------------------------ CLUSTERING
def calc_link(path_rmsd, path_link, link_method):
    print(f">>> Preparing 'clustering/{path_link.name}'...")
    rmsd_mat = np.load(path_rmsd)
    Z = linkage(rmsd_mat, link_method)
    np.save(path_link, Z, allow_pickle = False)
    print("...>>> Done.")

def calc_cluster(Z, t, label_criterion):
    return fcluster(Z, t = t, criterion = label_criterion)


# ------------------------------------------------------------------------------ CMAP
def calc_cmap(coords_frame):
    return distances.distance_array(coords_frame, coords_frame)

def calc_cmap_AS(coords0, coords1):
    return distances.distance_array(coords0, coords1)


# ------------------------------------------------------------------------------ PCA
def calc_pca(traj, pca_space_dir, pca_cumvar_dir, path_pca_pcomps):
    print(f">>> Preparing 'pca/{pca_space_dir.name}'...")

    print("...>>> Running PCA...")
    rec_pca = PCA(traj, select = "name CA")
    rec_pca.run()
    print("...>>> Done.")

    space = rec_pca.transform(traj.select_atoms("name CA"), 3)
    cumvar = rec_pca.cumulated_variance
    pcomps = rec_pca.p_components

    np.save(pca_space_dir, space, allow_pickle = False)
    np.save(pca_cumvar_dir, cumvar, allow_pickle = False)
    np.save(path_pca_pcomps, pcomps, allow_pickle = False)
    print(f"...>>> Saved to '{pca_space_dir.name}'.")


# ------------------------------------------------------------------------------ RAMA
def calc_rama(traj, path_rama):
    print(f">>> Preparing 'rama/{path_rama.name}'...")

    prot = traj.select_atoms("protein")
    rama = Ramachandran(prot).run()
    np.save(path_rama, rama.angles, allow_pickle = False)
    print("...>>> Done.")


# ------------------------------------------------------------------------------ RGYR
def calc_rgyr(traj, path_rgyr):
    print(f">>> Preparing 'rgyr/{path_rgyr.name}'...")

    rgyr = [traj.select_atoms("protein and name CA").radius_of_gyration() for _ in traj.trajectory[:]]
    np.save(path_rgyr, rgyr, allow_pickle = False)
    print("...>>> Done.")


# ------------------------------------------------------------------------------ RMSD
def calc_rmsd(coords, path_rmsd):
    print(f">>> Preparing 'rmsd/{path_rmsd.name}'...")
    n_frames = coords.shape[0]

    print(f"...>>> Calculating (non-redundant) all-to-all RMSD for {n_frames} frames...")
    result = [[rmsd(coords[i], coords[j]) for j in range(i)] for i in range(1, n_frames)]

    print(f"...>>> Fitting into a {n_frames}x{n_frames} matrix...")
    rmsd_map = np.zeros((n_frames, n_frames))
    for i,res in enumerate(result, 1):
        rmsd_map[i,:i] = rmsd_map[:i,i] = res

    np.save(path_rmsd, rmsd_map, allow_pickle = False)
    print(f"...>>> Saved to '{path_rmsd.name}'.")


# ------------------------------------------------------------------------------ RMSF
def calc_rmsf(traj, path_rmsf):
    print(f">>> Preparing 'rmsf/{path_rmsf.name}'...")

    ca_atoms = traj.select_atoms("protein and name CA")
    print("...>>> Calculating RMSF...")
    rmsf_CA = RMSF(ca_atoms).run()
    np.save(path_rmsf, rmsf_CA.rmsf, allow_pickle = False)
    print("...>>> Done.")


# //////////////////////////////////////////////////////////////////////////////
