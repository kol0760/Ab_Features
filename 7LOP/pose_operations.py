import pyrosetta.rosetta.core.pack.task as pyrosetta_task
import pyrosetta.rosetta.core.select as pyrosetta_select
from pyrosetta.rosetta.protocols import minimization_packing
from pyrosetta.rosetta.core import kinematics
import pyrosetta

def MutateResidue(pose, posi, amino, repack_radius=5.0, min_tolerance=0.01, max_iterations=200):
    """
    Function to pack rotamers and minimize around a specific residue in a given pose.
    
    Args:
        pose (Pose): The protein structure to be modified.
        posi (int): Residue index to mutate (get by "pose.pdb_info().pdb2pose(chain_id, res=int(feature['res']))").
        amino (str): Amino acid code for the residue to be designed (e.g. "A").
        scorefxn (ScoreFunction): The score function used for packing and minimization.
        repack_radius (float, optional): Radius for selecting neighbors for repacking. Default is 5.0 Å.
        min_tolerance (float, optional): Tolerance for the minimizer. Default is 0.01.
        max_iterations (int, optional): Maximum number of minimization iterations. Default is 200.
    """
    mut_posi = pyrosetta_select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)

    # Create a neighborhood selector around the mutation position
    nbr_selector = pyrosetta_select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)
    nbr_selector.set_distance(repack_radius)  
    not_design = pyrosetta_select.residue_selector.NotResidueSelector(mut_posi)

    # Task Factory setup
    tf = pyrosetta_task.TaskFactory()
    tf.push_back(pyrosetta_task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta_task.operation.IncludeCurrent())
    tf.push_back(pyrosetta_task.operation.NoRepackDisulfides())

    # Prevent repacking of the neighborhood residues
    prevent_repacking_rlt = pyrosetta_task.operation.PreventRepackingRLT()
    prevent_subset_repacking = pyrosetta_task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True)
    tf.push_back(prevent_subset_repacking)

    # Restrict non-design residues to repacking only
    tf.push_back(pyrosetta_task.operation.OperateOnResidueSubset(
        pyrosetta_task.operation.RestrictToRepackingRLT(), not_design))
    
    # Restrict the mutation position to the specified amino acid
    aa_to_design = pyrosetta_task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(pyrosetta_task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

    # Pack Rotamers Mover
    packer = minimization_packing.PackRotamersMover()
    packer.task_factory(tf)
    packer.apply(pose)


    # Set up MoveMap for minimization
    movemap = kinematics.MoveMap()
    movemap.set_bb(False)  
    movemap.set_chi(False) 
    for res_index in pyrosetta_select.get_residues_from_subset(nbr_selector.apply(pose)):
        movemap.set_bb(res_index, True)
        movemap.set_chi(res_index, True)
    scorefxn = pyrosetta.create_score_function('ref2015')    
    # Minimization Mover
    min_mover = minimization_packing.MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)
    min_mover.min_type("lbfgs_armijo_nonmonotone") 
    min_mover.tolerance(min_tolerance)
    min_mover.max_iter(max_iterations)
    min_mover.apply(pose)
