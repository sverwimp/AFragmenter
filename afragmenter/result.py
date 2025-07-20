from pathlib import Path
from typing import Optional, Union, Tuple

from matplotlib import image, axes

from .sequence_reader import SequenceReader
from afragmenter import format_result, plotting
from .structure_viewer import view_py3Dmol


FilePath = Union[str, Path]


class ClusteringResult:
    """
    Class to hold the clustering results of the AFragmenter.

    Parameters:
    - pae_matrix (np.ndarray): The Predicted Aligned Error matrix.
    - cluster_intervals (dict): The cluster intervals obtained from the clustering results. Each key represents a cluster number,
                                and the value is a list of tuples containing the start and end indices of the cluster intervals.
    - params (dict): The parameters used for clustering.
    - sequence_reader (Optional[SequenceReader]): The SequenceReader object used to read the sequence file. This is automatically created
                                                 when the print_fasta or save_fasta methods are called
    """
    def __init__(self, pae_matrix, cluster_intervals, params, sequence_reader=None):
        self.pae_matrix = pae_matrix
        self.cluster_intervals = cluster_intervals
        self.params = params
        self.sequence_reader = sequence_reader

    def plot_pae(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        """
        Plot the Predicted Aligned Error matrix as a heatmap.

        Parameters:
        - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

        Returns:
        - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
        """
        kwargs.setdefault("colorbar_label", "Predicted Aligned Error (Å)")
        image, ax = plotting.plot_matrix(self.pae_matrix, **kwargs)
        ax.set_xlabel("Scored residue")
        ax.set_ylabel("Aligned residue")

        return image, ax

    def plot_pae(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        """
        Plot the Predicted Aligned Error matrix as a heatmap.

        Parameters:
        - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

        Returns:
        - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
        """
        kwargs.setdefault("colorbar_label", "Predicted Aligned Error (Å)")
        image, ax = plotting.plot_matrix(self.pae_matrix, **kwargs)
        ax.set_xlabel("Scored residue")
        ax.set_ylabel("Aligned residue")

        return image, ax

    def plot_result(self, **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
        """
        Plot the clustering results on top of the Predicted Aligned Error matrix.

        Parameters:
        - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

        Returns:
        - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
        """
        
        kwargs.setdefault("colorbar_label", "Predicted Aligned Error (Å)")
        image, ax = plotting.plot_cluster_intervals(self.pae_matrix, self.cluster_intervals, **kwargs)
        ax.set_xlabel("Scored residue")
        ax.set_ylabel("Aligned residue")

        return image, ax

    def print_result(self, base_name: Optional[str] = None, format: str = 'auto', delimiter: str = ',') -> None:
        """
        Print the clustering results in a table format. Will use rich if the output is a terminal or a jupyter notebook, else will use csv.

        Parameters:
        - base_name (str, optional): the base name to use for the cluster intervals.
        - format (str, optional): The format to use. [auto, csv, rich]
        - delimiter (str, optional): The delimiter to use for the csv format.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        
        # not 'if not base_name' because base_name can be an empty string such that only the cluster number is used as name
        if base_name is None:
            base_name = self.sequence_reader.name if self.sequence_reader else None

        format_result.print_result(self.cluster_intervals, format=format, delimiter=delimiter, base_name=base_name)

    
    def save_result(self, output_file: FilePath, base_name: Optional[str] = None, format: str = 'csv', delimiter: str = ',') -> None:
        """
        Save the clustering results in a table format to a file.

        Parameters:
        - output_file (FilePath): The path to save the output file.
        - base_name (str, optional): the base name to use for the cluster intervals.
        - format (str, optional): The format to use. [auto, csv, rich]
        - delimiter (str, optional): The delimiter to use for the csv format.

        Returns:
        - None

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        
        if base_name is None:
            base_name = self.sequence_reader.name if self.sequence_reader else None

        format_result.save_result(self.cluster_intervals, output_file, base_name=base_name, format=format, delimiter=delimiter)
    
    
    def py3Dmol(self, structure_file: str,
                color_palette: Optional[list] = None,
                style: str = 'cartoon',
                displace_domains: bool = False,
                tether_strength: float = 0.1,
                repulse_strength: float = 100,
                steps: int = 1000,
                dt: float = 0.01,
                **kwargs) -> 'py3Dmol.view':
        """
        Visualize the 3D structure of the protein using py3Dmol. Color the residues based on the clusters.

        Parameters:
        - structure_file (str): The path to the PDB or mmcif file.
        - color_palette (list, optional): A list of colors to use for the clusters, expects color names or hex-codes.
        - style (str, optional): The style to use. Defaults to 'cartoon'.
        - displace_domains (bool, optional): Whether to displace the domains based on the clusters. Defaults to False.
        - tether_strength (float, optional): The spring constant pulling each group toward its original center. Defaults 
                                             to 0.1. Only used if displace_domains is True.
        - repulse_strength (float, optional): The constant for the inverse-square repulsion between groups. Defaults to 100.
                                              Only used if displace_domains is True.
        - steps (int, optional): The number of simulation iterations. Defaults to 1000. Only used if displace_domains is True.
        - dt (float, optional): The time step for the simulation. Defaults to 0.01. Only used if displace_domains is True.
        - **kwargs: Additional keyword arguments to be passed to the py3Dmol.view function.

        Returns:
        - py3Dmol.view: The py3Dmol viewer object.

        Raises:
        - ImportError: If the py3Dmol library is not installed.
        - ValueError: If the cluster intervals are not defined.
        - ValueError: If the structure file is not a PDB or mmCIF file.
        """
        
        if color_palette is None:
            from .structure_viewer import COLOR_PALETTE
            color_palette = COLOR_PALETTE

        view = view_py3Dmol(structure_file,
                            self.cluster_intervals,
                            displace_domains=displace_domains,
                            color_palette=color_palette,
                            style=style,
                            tether_strength=tether_strength,
                            repulse_strength=repulse_strength,
                            steps=steps,
                            dt=dt,
                            **kwargs)
        return view


    def _format_fasta_sequences(self, parsed_sequences: dict, header_name: str, width: int) -> str: # type: ignore (generator)
        """
        Generator function to format sequences in FASTA format.

        Parameters:
        - parsed_sequences (dict): The parsed sequences.
        - prefix (str): The name to add to the sequence headers.
        - width (int): The width of the sequence lines.

        Yields:
        - str: The formatted FASTA sequence.
        """
        for i, seq in parsed_sequences.items():
            interval = self.cluster_intervals[i]
            interval_str = "_".join([f"{start+1}-{end+1}" for start, end in interval])
            name = f"{header_name}_{i+1}" if header_name else f"{i+1}"
            yield f">{name} {interval_str}"
            for j in range(0, len(seq), width):
                yield seq[j:j+width]


    def _get_parsed_sequence(self, sequence_file: FilePath) -> SequenceReader:
        """
        Get the parsed sequences from the sequence file.
        
        Parameters:
        - sequence_file (FilePath): The path to the sequence file.

        Returns:
        - dict: A dictionary where the keys are the cluster indices and the values are the corresponding protein sequences.
        """

        if self.sequence_reader is None:
            self.sequence_reader = SequenceReader(sequence_file, seq_length=self.pae_matrix.shape[0])
        parsed_sequences = self.sequence_reader.clusters_to_fasta(self.cluster_intervals)

        return parsed_sequences


    def print_fasta(self, sequence_file: FilePath, header_name: Optional[str] = None, width: int = 60) -> None:
        """
        Print the sequences corresponding to each cluster in FASTA format.

        Parameters:
        - sequence_file (FilePath): The path to the sequence file.
        - header_name (str, optional): The name to add to the sequence headers. Defaults to self.sequence_reader.name.
        - width (int, optional): The width of the sequence lines. Defaults to 60.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        parsed_sequences = self._get_parsed_sequence(sequence_file)
        
        if header_name is None:
            header_name = self.sequence_reader.name

        for line in self._format_fasta_sequences(parsed_sequences, header_name, width):
            print(line)


    def save_fasta(self, sequence_file: FilePath, output_file: FilePath, header_name: Optional[str] = None, width: int = 60) -> None:
        """
        Save the sequences corresponding to each cluster in FASTA format to a file.

        Parameters:
        - sequence_file (FilePath): The path to the sequence file.
        - output_file (FilePath): The path to save the output file.
        - header_name (str, optional): The name to add to the sequence headers. Defaults to self.sequence_reader.name.
        - width (int, optional): The width of the sequence lines. Defaults to 60.

        Raises:
        - ValueError: If the cluster intervals are not defined.
        """
        parsed_sequences = self._get_parsed_sequence(sequence_file)
        
        if header_name is None:
            header_name = self.sequence_reader.name

        output_file = Path(output_file)

        with output_file.open('w') as f:
            for line in self._format_fasta_sequences(parsed_sequences, header_name, width):
                f.write(line + "\n")

