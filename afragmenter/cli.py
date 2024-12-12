import sys
import importlib.metadata
from pathlib import Path

import rich_click as click

from .afragmenter import AFragmenter


click.rich_click.OPTION_GROUPS = {
    "afragmenter": [
        {
            "name": "Input",
            "options": [
                "--structure_file",
                "--json",
            ]
        },
        {
            "name": "Clustering options",
            "options": [
                "--resolution",
                "--objective_function",
                "--n_iterations",
            ]
        },
        {
            "name": "Fine tuning",
            "options": [
                "--threshold",
                "--min_size",
            ]
        },
        {
            "name": "Outputs",
            "options": [
                "--plot_results",
                "--output_fasta",
            ]
        }
    ],
}


@click.command()
@click.help_option("--help", "-h")
@click.version_option(importlib.metadata.version("afragmenter"), "--version", "-V")
@click.option("--structure_file", 
              "-s", type=click.Path(exists=True), 
              required=False, 
              help="Path to a PDB or mmCIF file containing the protein structure"
)
@click.option("--json", 
              "-j" , 
              type=click.Path(exists=True), 
              required=True, 
              help="Path to the AlphaFold json file containing the PAE data"
)
@click.option("--resolution", 
              "-r", 
              type=click.FloatRange(min=0.0, min_open=True), 
              help="Resolution used with Leiden clustering (default = 0.8 for modularity, and 0.3 for CPM)"
)
@click.option("--objective_function", 
              "-f", 
              type=click.Choice(["modularity", "CPM"], case_sensitive=False), 
              default="modularity",
              show_default=True,
              help="Objective function for Leiden clustering (not case sensitive)"
)
@click.option("--n_iterations", 
              "-n", 
              type=click.IntRange(),
              default=1_000, 
              show_default=True,
              help="Number of iterations for Leiden clustering. Negative values will run the algorithm until a stable iteration is reached",
)
@click.option("--threshold",
              type=click.FloatRange(min=0.0, max=31.75), 
              default=5, 
              show_default=True,
              help="Threshold for the sigmoid function used to transform the PAE values into graph edge weights"
)
@click.option("--min_size",
              type=click.IntRange(min=0), 
              default=0, 
              show_default=True,
              help="Minimum cluster size. Maximum size is equal to the number of residues in the protein"
)
@click.option("--plot_results",
              type=click.Path(path_type=Path, dir_okay=False, writable=True),
              default=None,
              help="Path to save the results plot"
)
@click.option("--output_fasta",
              type=click.Path(path_type=Path, dir_okay=False, writable=True),
              default=None,
              help="Path to save the output fasta file (requires --structure_file)"
)
def main(structure_file: Path, 
         json: Path,
         resolution: float,
         objective_function: str,
         n_iterations: int,
         threshold: float,
         min_size: int,
         plot_results: Path,
         output_fasta: Path):
    
    afragmenter = AFragmenter(pae_matrix=json, 
                              threshold=threshold)

    afragmenter.cluster(resolution=resolution, 
                        n_iterations=n_iterations, 
                        objective_function=objective_function, 
                        min_size=min_size)

    afragmenter.print_results()

    if plot_results:
        _, ax = afragmenter.plot_results()
        fig = ax.get_figure()
        fig.savefig(plot_results)

    if output_fasta:
        if not structure_file:
            raise click.BadOptionUsage("output_fasta", "The --structure_file option is required when using --output_fasta")
        afragmenter.save_fasta(sequence_file=structure_file, output_file=output_fasta)


if __name__ == "__main__":
    main()
    sys.exit(0)