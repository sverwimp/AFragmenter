import sys
import importlib.metadata
from pathlib import Path

import rich_click as click
from rich_click import Option

from .afragmenter import AFragmenter
from .afdb_client import fetch_afdb_data


click.rich_click.USE_RICH_MARKUP = True
click.rich_click.OPTION_GROUPS = {
    "afragmenter": [
        {
            "name": "Input",
            "options": [
                "--json",
                "--structure",
                "--afdb",
            ]
        },
        {
            "name": "Clustering options",
            "options": [
                "--resolution",
                "--objective-function",
                "--n-iterations",
            ]
        },
        {
            "name": "Fine tuning",
            "options": [
                "--threshold",
                "--min-size",
                "--no-merge",
            ]
        },
        {
            "name": "Outputs",
            "options": [
                "--plot-results",
                "--output-fasta",
            ]
        }
    ],
}


class EitherRequired(Option):
    def __init__(self, *args, **kwargs):
        self.either_required = kwargs.pop("either_required", [])
        super().__init__(*args, **kwargs)
    
    def handle_parse_result(self, ctx, opts, args):
        if not hasattr(ctx, 'either_required_message_printed'):
            ctx.either_required_message_printed = False
        
        # Print a message if both --json and --afdb are provided, once
        if set(self.either_required).intersection(opts) and self.name in opts:
            if not ctx.either_required_message_printed:
                print("Both --json and --afdb were provided. Using --afdb")
                ctx.either_required_message_printed = True
        
        if not any(opt in opts for opt in self.either_required + [self.name]):
            raise click.UsageError(f"Either --{self.name} or --{' or --'.join(self.either_required)} is required.")
            
        return super().handle_parse_result(ctx, opts, args)


@click.command()
@click.help_option("--help", "-h")
@click.version_option(importlib.metadata.version("afragmenter"), "--version", "-v")
@click.option("--structure", 
              "-s", type=click.Path(exists=True, dir_okay=False, readable=True), 
              required=False, 
              help="Path to a PDB or mmCIF file containing the protein structure and sequence"
)
@click.option("--json", 
              "-j" , 
              type=click.Path(exists=True, dir_okay=False, readable=True), 
              cls=EitherRequired,
              either_required=["afdb"],
              help="Path to the AlphaFold json file containing the PAE data [red]\[required: either --json or --afdb][/]"
)
@click.option("--afdb",
              type=str,
              cls=EitherRequired,
              either_required=["json"],
              help="Uniprot identifier to fetch data from the AlphaFold database [red]\[required: either --json or --afdb][/]"
)
@click.option("--resolution", 
              "-r", 
              type=click.FloatRange(min=0.0, min_open=True), 
              help="Resolution used with Leiden clustering [dim]\[default: 0.8 for modularity, 0.3 for CPM][/]",
)
@click.option("--objective-function", 
              "-f", 
              type=click.Choice(["modularity", "CPM"], case_sensitive=False), 
              default="modularity",
              show_default=True,
              help="Objective function for Leiden clustering (not case sensitive)"
)
@click.option("--n-iterations", 
              "-n",
              type=click.INT,
              default=10_000, 
              show_default=True,
              help="Number of iterations for Leiden clustering. Negative values will run the algorithm until a stable iteration is reached",
)
@click.option("--threshold",
              "-t",
              type=click.FloatRange(min=0.0, max=31.75), 
              default=5, 
              show_default=True,
              help="Contrast thresold for the PAE values. This is a soft cut-off point to increase the contrast between low and high PAE values. Values near the threshold will transition more smoothly."
              #"Threshold for the sigmoid function used to transform the PAE values into graph edge weights"
)
@click.option("--min-size",
              "-m",
              type=click.IntRange(min=0), 
              default=10, 
              show_default=True,
              help="Minimum size of partition to be considered. Attempts to merge partitions that are too small with adjecent larger ones. "
                   "Set to 0 to keep all partitions."
              #"Minimum cluster size. Maximum size is equal to the number of residues in the protein"
)
@click.option("--no-merge",
              is_flag=True,
              default=False,
              help="Do not attempt to merge small partitions with larger paritions, just discard them",
)
@click.option("--plot-results",
              type=click.Path(path_type=Path, dir_okay=False, writable=True),
              default=None,
              help="Path to save the results plot"
)
@click.option("--output-fasta",
              type=click.Path(path_type=Path, dir_okay=False, writable=True),
              default=None,
              help="Path to save the output fasta file (requires --structure-file)"
)
def afragmenter(structure: Path, 
         json: Path,
         afdb: str,
         resolution: float,
         objective_function: str,
         n_iterations: int,
         threshold: float,
         min_size: int,
         no_merge: bool,
         plot_results: Path,
         output_fasta: Path):
    
    if afdb:
        json, structure = fetch_afdb_data(afdb)
    
    afragmenter = AFragmenter(pae_matrix=json, 
                              threshold=threshold)

    attempt_merge = not no_merge
    afragmenter.cluster(resolution=resolution, 
                        n_iterations=n_iterations, 
                        objective_function=objective_function, 
                        min_size=min_size,
                        attempt_merge=attempt_merge)

    afragmenter.print_result()

    if plot_results:
        _, ax = afragmenter.plot_result()
        fig = ax.get_figure()
        fig.savefig(plot_results)

    if output_fasta:
        if not structure:
            raise click.BadOptionUsage("output_fasta", "The --structure-file option is required when using --output-fasta")
        afragmenter.save_fasta(sequence_file=structure, output_file=output_fasta)


if __name__ == "__main__":
    afragmenter()
    sys.exit(0)