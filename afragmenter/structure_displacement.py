from io import StringIO
from typing import Union

import numpy as np
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser


def _compute_group_center(group: list) -> np.ndarray:
    """
    Compute the center of a group of residues using the coordinates of the CA atoms.
    If a residue does not have a CA atom, it is skipped.

    Parameters:
    - group: A list of Bio.PDB.Residue.Residue objects.

    Returns:
    - The center of the group as a numpy array, or None if no CA atoms are found.
    """
    coords = []
    for residue in group:
        if "CA" in residue:
            coords.append(residue["CA"].get_coord())
    if not coords:
        return None
    return np.mean(coords, axis=0)


def _simulate_with_tether(initial_centers: list, 
                          tether_strength: float = 0.1, 
                          repulse_strength: float = 100, 
                          steps: int = 1000, 
                          dt: float = 0.01) -> np.ndarray:
    """
    Given a list of original centers (as np.array), simulate their displacement using:
      - A tether force that pulls each center toward its original position.
      - A repulsive force that pushes centers apart.
    Returns the new positions for each center.

    Parameters:
    - initial_centers (list[np.ndarray]): List of initial center positions as np.array.
    - tether_strength (float): Spring constant for the tether force.
    - repulse_strength (float): Constant for the inverse-square repulsion between centers.
    - steps (int): Number of simulation steps.
    - dt (float): Time step for the simulation.

    Returns:
    - A numpy array of the new positions for each center.

    Note: This function uses Euler integration for simplicity.
    """
    positions = np.array(initial_centers, dtype=float)
    velocities = np.zeros_like(positions)
    
    for _ in range(steps):
        forces = np.zeros_like(positions)
        
        # Tether force: pull each center back to its original position.
        for i in range(len(positions)):
            forces[i] += -tether_strength * (positions[i] - initial_centers[i])
        
        # Repulsive forces: push centers away from each other.
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                diff = positions[i] - positions[j]
                distance = np.linalg.norm(diff)
                if distance < 1e-5:
                    distance = 1e-5  # avoid division by zero
                force_magnitude = repulse_strength / (distance ** 2)
                force_vector = (diff / distance) * force_magnitude
                forces[i] += force_vector
                forces[j] -= force_vector
        
        # Euler integration update.
        velocities += forces * dt
        positions += velocities * dt
        
    return positions


def displace_structure(input_content: str, 
                       moving_group_intervals: Union[list, dict] = None,
                       file_format: str = 'pdb', 
                       tether_strength: float = 0.1, 
                       repulse_strength: float = 100, 
                       steps: int = 1000, 
                       dt: float = 0.01):
    """
    Reads a structure from input_file (PDB or mmCIF), displaces specified moving groups
    using tether (keeping them near their original centers) and repulsive forces, and
    returns the modified structure as a string.

    Parameters:
    - input_content (str): A string containing the content of the structure file.
    - moving_group_intervals (list): List where each element represents one moving group. Each group is defined
                                     as a list of one or more tuples (start, end) that specify intervals of residue
                                     indices (based on the ordered list of standard residues) to be moved.
                                     A dictionary can also be passed, in which case the values are used as the
                                     moving group intervals.
    - file_format (str): Either 'pdb' or 'mmcif'.
    - tether_strength (float): Spring constant pulling each group toward its original center.
    - repulse_strength (float): Constant for the inverse-square repulsion between groups.
    - steps (int): Number of simulation iterations.
    - dt (float): Time step for the simulation.

    Returns:
    - A string containing the modified structure in the same format as the input.
    
    Example usage:
    ```
    input_content = open("example.pdb").read()
    moving_group_intervals = [
        [(5, 15), (30, 40)],  # Group 0 (discontinuous intervals)
        [(15, 30)]            # Group 1 (continuous interval)
    ]
    modified_structure = displace_structure(input_file, moving_group_intervals)
    ```
    """
    # Choose the correct parser.
    if file_format.lower() in ['mmcif', 'cif']:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    # Parse the input content.
    file_obj = StringIO(input_content)
    structure = parser.get_structure("structure", file_obj)
    
    # Extract all standard residues (assumes one chain; modify if needed).
    all_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == " ":  # Only standard residues
                    all_residues.append(residue)
    
    # Create mapping for residues in moving groups.
    # 'moving_groups' will be a list where each element is a list of residues (which may be discontinuous).
    # 'moving_map' will map a residue's id to its group index.
    moving_groups = []
    moving_map = {}

    if isinstance(moving_group_intervals, dict):
        moving_group_intervals = list(moving_group_intervals.values())
    
    if moving_group_intervals is None:
        # If no groups specified, treat all residues as one moving group.
        moving_groups.append(all_residues)
        for residue in all_residues:
            moving_map[id(residue)] = 0
    else:
        for group_idx, intervals in enumerate(moving_group_intervals):
            # Allow a single tuple to be passed instead of a list.
            if isinstance(intervals, tuple):
                intervals = [intervals]
            group = []
            for (start, end) in intervals:
                # Append residues in the interval [start, end) from the ordered list.
                for residue in all_residues[start:end]:
                    group.append(residue)
                    moving_map[id(residue)] = group_idx
            moving_groups.append(group)
    
    # Compute the original centers for each moving group.
    original_centers = []
    for group in moving_groups:
        center = _compute_group_center(group)
        if center is None:
            raise ValueError("A moving group has no CA atoms to compute its center.")
        original_centers.append(center)
    
    # Run the simulation to obtain new centers.
    new_centers = _simulate_with_tether(original_centers,
                                       tether_strength=tether_strength,
                                       repulse_strength=repulse_strength,
                                       steps=steps, 
                                       dt=dt)
    
    # Compute displacement vectors for each moving group.
    displacements = [new_centers[i] - original_centers[i] for i in range(len(moving_groups))]
    
    # Apply the displacement to residues that are in moving groups.
    for residue in all_residues:
        if id(residue) in moving_map:
            group_idx = moving_map[id(residue)]
            disp = displacements[group_idx]
            for atom in residue:
                atom.set_coord(atom.get_coord() + disp)
        # Residues not in any moving group remain static.
    
    # Write the modified structure to a string.
    output_buffer = StringIO()
    if file_format.lower() in ['mmcif', 'cif']:
        writer = MMCIFIO()
    else:
        writer = PDBIO()
    writer.set_structure(structure)
    writer.save(output_buffer)
    modified_structure_str = output_buffer.getvalue()
    
    return modified_structure_str
