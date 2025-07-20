import os
from typing  import TYPE_CHECKING

from .sequence_reader import SequenceReader
from .structure_displacement import displace_structure

COLOR_PALETTE = ['red', 'blue', 'green', 'yellow', 'purple', 'orange', 'cyan', 'magenta', 
                'lime', 'pink', 'teal', 'lavender', 'brown', 'apricot', 'maroon', 'mint', 'olive', 
                'beige', 'navy', 'grey', 'white', 'black']

if TYPE_CHECKING:
    import py3Dmol  # type: ignore # This import is only for type checking and should not be executed at runtime.


def color_view_by_domain(view: 'py3Dmol.view',  # type: ignore
                    cluster_intervals: dict, 
                    color_palette: list = COLOR_PALETTE, 
                    style: str = 'cartoon') -> None:
    """
    Color the structure in the view according to the domain clusters provided in cluster_intervals.
    Each cluster will be colored with a different color from the color_range.

    Parameters:
    - view (py3Dmol.view): The view to color.
    - cluster_intervals (dict): Dictionary where keys are cluster indices and values are lists of (start, end) residue indices.
    - color_palette (list): List of colors to use for the clusters.
    - style (str): Style to apply.
    """
    if len(cluster_intervals) > len(color_palette):
        print("Warning: More clusters than available colors. Some clusters will have the same color.")
    
    for cluster, ranges in cluster_intervals.items():
        col = color_palette[cluster % len(color_palette)]
        for start, end in ranges:
            view.setStyle({'resi': f'{start+1}-{end+1}'}, {style: {'color': col}})
    


def view_py3Dmol(structure_file: str, 
                 cluster_intervals: dict, 
                 displace_domains: bool = False, 
                 color_palette: list = COLOR_PALETTE,
                 style: str = 'cartoon',
                 **kwargs) -> 'py3Dmol.view': # type: ignore
    """
    Create a 3D view of the structure in structure_file, coloring the domains according to cluster_intervals.

    Parameters:
    - structure_file (str): Path to the structure file (PDB or mmCIF) or the structure content.
    - cluster_intervals (dict): Dictionary where keys are cluster indices and values are lists of (start, end) residue indices.
    - displace_domains (bool): Whether to displace the domains using tether and repulsive forces.
    - color_palette (list): List of colors to use for the clusters.
    - style (str): Style to apply.
    
    Returns:
    - The py3Dmol.view object.
    """
    try:
        # Check if py3Dmol is installed for view_py3Dmol function
        import py3Dmol  # type: ignore
    except ImportError:
        raise ImportError(
            "The py3Dmol library is required for the py3Dmol function. "
            "Please install it using 'pip install py3Dmol'."
            )

    # Read the structure and determine the file format.
    if os.path.isfile(structure_file):
        content = open(structure_file).read()
    else:
        content = structure_file
    file_format = SequenceReader.determine_file_format(content)

    # Displace the domains if requested.
    tether_strength = kwargs.pop('tether_strength', 0.1)
    repulse_strength = kwargs.pop('repulse_strength', 100)
    steps = kwargs.pop('steps', 1000)
    dt = kwargs.pop('dt', 0.01)
    if displace_domains:
        content = displace_structure(content, 
                                     moving_group_intervals=cluster_intervals, 
                                     file_format=file_format,
                                     tether_strength=tether_strength,
                                     repulse_strength=repulse_strength,
                                     steps=steps,
                                     dt=dt)

    # Create the view.
    kwargs.setdefault('width', 800)
    kwargs.setdefault('height', 600)
    view = py3Dmol.view(**kwargs)
    if file_format == 'mmcif':
        view.addModel(content, 'mmcif')
    elif file_format == 'pdb':
        view.addModel(content, 'pdb')
    else:
        raise ValueError(f'Unsupported file format: {file_format}, please provide a PDB or mmCIF file.')
    
    view.setStyle({style: {'color': 'grey'}}) # Default color.
    # Color by domain.
    color_view_by_domain(view, cluster_intervals, color_palette, style)

    return view
