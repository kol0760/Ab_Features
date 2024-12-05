from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.core.conformation import disulfide_bonds
from pyrosetta.rosetta.utility import vector1_std_pair_unsigned_long_unsigned_long_t
from tqdm import tqdm
from scipy.spatial.distance import cdist

import numpy as np

def hbond_analysis(pose, num_residues, interaction_matrix):
    # H-bond interaction analysis
    hbond_set = pose.get_hbonds()
    hbond_threshold = 4
    for hbond in range(1, hbond_set.nhbonds() + 1):
        hbond_obj = hbond_set.hbond(hbond)
        donor_res = hbond_obj.don_res()
        acceptor_res = hbond_obj.acc_res()
        H_bond_distance = hbond_obj.get_HAdist(pose)
        don_is_backbone = hbond_obj.don_hatm_is_backbone()
        acc_is_backbone = hbond_obj.acc_atm_is_backbone()
        # H_bond_ = round((hbond_threshold - H_bond_distance), 3)
        H_bond_ = 1
        H_bond_interaction = interaction_matrix[donor_res - 1][acceptor_res - 1]
        symmetric_interaction = interaction_matrix[acceptor_res - 1][donor_res - 1]

        if don_is_backbone and acc_is_backbone:
            H_bond_interaction[2] += H_bond_
            symmetric_interaction[2] += H_bond_
        elif (not don_is_backbone) and (not acc_is_backbone):
            H_bond_interaction[0] += H_bond_
            symmetric_interaction[0] += H_bond_
        else:
            H_bond_interaction[1] += H_bond_
            symmetric_interaction[1] += H_bond_

def vdw_analysis(pose, num_residues, interaction_matrix):
    # van der Waals interaction analysis
    vdw_threshold = 5
    all_coords = []
    atom_indices = []  #  (residue index, atom index) 
    for i in range(num_residues):
        res = pose.residue(i + 1)
        for atom in range(1, res.natoms() + 1):
            if res.atom_type(atom).element() == 'C':  
                all_coords.append(np.array(res.xyz(atom)))
                atom_indices.append((i, atom))  

    all_coords = np.array(all_coords)
    dist_matrix = cdist(all_coords, all_coords)
    mask = dist_matrix <= vdw_threshold

    for i in range(len(mask)):
        for j in range(i + 1, len(mask)): 
            if mask[i, j]:
                res1_idx, atom1_idx = atom_indices[i]
                res2_idx, atom2_idx = atom_indices[j]
                if res1_idx != res2_idx:
                    res1 = pose.residue(res1_idx + 1)
                    res2 = pose.residue(res2_idx + 1)

                    res1_is_backbone = res1.atom_is_backbone(atom1_idx)
                    res2_is_backbone = res2.atom_is_backbone(atom2_idx)

                    vdw_interaction = interaction_matrix[res1_idx][res2_idx]
                    symmetric_interaction = interaction_matrix[res2_idx][res1_idx]
                    # vdw_ = round((vdw_threshold - distance), 3)
                    vdw_ = 1

                    if res1_is_backbone and res2_is_backbone:
                        vdw_interaction[5] += vdw_
                        symmetric_interaction[5] += vdw_
                    elif (not res1_is_backbone) and (not res2_is_backbone):
                        vdw_interaction[3] += vdw_
                        symmetric_interaction[3] += vdw_
                    else:
                        vdw_interaction[4] += vdw_
                        symmetric_interaction[4] += vdw_

def disulfide_analysis(pose,num_residues, interaction_matrix):
    # Disulfide bond analysis

    disulfide_threshold = 2.5
    disulfides = vector1_std_pair_unsigned_long_unsigned_long_t()
    disulfide_bonds(pose.conformation(), disulfides)
    for disulfide_pair in disulfides:
        res1_idx = disulfide_pair[0]
        res2_idx = disulfide_pair[1]
        res1 = pose.residue(res1_idx)
        res2 = pose.residue(res2_idx)
        distance = res1.atom('SG').xyz().distance(res2.atom('SG').xyz())
        if distance < disulfide_threshold:
            # disulfide_ = round((disulfide_threshold - distance), 3)
            disulfide_ = 1
            interaction_matrix[res1_idx - 1][res2_idx - 1][7] += disulfide_
            interaction_matrix[res2_idx - 1][res1_idx - 1][7] += disulfide_

def pipi_analysis(pose, num_residues, interaction_matrix):
    # Pi-pi interaction analysis
    def pipi_interaction(center1, center2, normal1, normal2):
        distance = np.linalg.norm(center1 - center2)
        cos_theta = np.dot(normal1, normal2)
        angle = np.arccos(cos_theta) * (180 / np.pi)
        if angle > 90:
            angle = 180 - angle
        vector_between_centers = center2 - center1
        projection = center2 - np.dot(vector_between_centers, normal1) * normal1
        offset = np.linalg.norm(projection - center1)
        distance_ok = 3.5 <= distance <= 5.5
        offset_ok = offset <= 3
        if 30 <= angle <= 36 or 60 <= angle <= 70:
            offset_ok = offset < 1
        elif angle >= 70:
            offset_ok = offset <= 5
        angle_ok = angle <= 36 or 60 <= angle
        if distance_ok and angle_ok and offset_ok:
            return distance
        else:
            return 0

    def extract_ring_centers_and_normals(ring_residue):
        ring_atoms = ring_residue.type().ring_atoms()
        ring_centers_and_normals = []
        for ring_index in range(ring_atoms.l(), ring_atoms.u() + 1):
            atom_indices = ring_atoms[ring_index]
            coords = np.array([list(ring_residue.xyz(atom_index)) for atom_index in atom_indices])
            center = np.mean(coords, axis=0)
            v1 = coords[1] - coords[0]
            v2 = coords[2] - coords[0]
            normal_vector = np.cross(v1, v2)
            normal_vector /= np.linalg.norm(normal_vector)
            ring_centers_and_normals.append((center, normal_vector))
        return ring_centers_and_normals

    aromatic_residues = ["TYR", "PHE", "TRP"]
    pipi_threshold = 7
    Aromatic_residues = []
    for i in range(1, num_residues + 1):
        residue = pose.residue(i)
        if residue.name3() in aromatic_residues:
            Aromatic_residues.append(residue)

    residue_ring_data = []
    for residue in Aromatic_residues:
        centers_and_normals = extract_ring_centers_and_normals(residue)
        residue_ring_data.append((residue, centers_and_normals))

    for i, (residue1, rings1) in enumerate(residue_ring_data):
        for j, (residue2, rings2) in enumerate(residue_ring_data):
            if i >= j:
                continue
            for center1, normal1 in rings1:
                for center2, normal2 in rings2:
                    if pipi_interaction(center1, center2, normal1, normal2) != 0:
                        res1_idx = residue1.seqpos()
                        res2_idx = residue2.seqpos()
                        # pipi_ = round((pipi_threshold - pipi_interaction(center1, center2, normal1, normal2)), 3)
                        pipi_ = 1
                        interaction_matrix[res1_idx - 1][res2_idx - 1][6] += pipi_
                        interaction_matrix[res2_idx - 1][res1_idx - 1][6] += pipi_


def compute_sin_scores(pose):
    # Initialize PyRosetta
    # init("-mute core basic protocols")

    # Load PDB file
    # pose = pose_from_pdb(pdb_file)
    num_residues = pose.total_residue()

    # Initialize interaction matrix
    interaction_matrix = [[[0 for _ in range(8)] for _ in range(num_residues)] for _ in range(num_residues)]

    # Perform interaction analyses
    # pbar.set_description("Analyzing H-bond interactions")
    hbond_analysis(pose, num_residues, interaction_matrix)
    # pbar.update(1)

    # pbar.set_description("Analyzing van der Waals interactions")
    vdw_analysis(pose, num_residues, interaction_matrix)
    # pbar.update(1)

    # pbar.set_description("Analyzing Disulfide interactions")
    disulfide_analysis(pose, num_residues, interaction_matrix)
    # pbar.update(1)

    # pbar.set_description("Analyzing Pi-Pi interactions")
    pipi_analysis(pose, num_residues, interaction_matrix)
    # pbar.update(1)

    # Calculate SIN scores
    SIN_weights = np.array([0.20, 0.10, 0.05, 0.10, 0.05, 0.025, 0.50, 0.50])/1.525
    interaction_matrix_sum = np.tensordot(interaction_matrix, SIN_weights, axes=([2], [0]))
    interaction_matrix_sum = np.round(interaction_matrix_sum, 3)

    # SIN_scores = interaction_matrix_sum.sum(axis=1)
    # max_score = np.max(SIN_scores)
    # min_score = np.min(SIN_scores)
    # SIN_scores_normalized = np.round((SIN_scores - min_score) / (max_score - min_score), 6)

    # return SIN_scores_normalized,interaction_matrix_sum, (max_score - min_score)
    return interaction_matrix_sum


# pdb_file = "ref.pdb"
# SIN_scores_normalized, score_range = compute_sin_scores(pdb_file)
# print("Normalized SIN scores:", SIN_scores_normalized)
# print("Score range:", score_range)
