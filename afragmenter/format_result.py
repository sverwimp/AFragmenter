
from io import StringIO
import csv
from typing import Union
from pathlib import Path
import sys

from rich.table import Table
from rich.console import Console


FilePath = Union[str, Path]


def format_csv(intervals: dict, delimiter: str = ',') -> str:
    """
    Format the cluster_interval results as a CSV string.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - delimiter (str): The delimiter to use. Default is ','.

    Returns:
    - str: The CSV formatted output.
    """
    output = StringIO()
    writer = csv.writer(output, delimiter=delimiter)
    writer.writerow(["domain", "nres", "chopping"])

    for i, cluster in enumerate(intervals.values()):
        chopping = '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])
        nres = sum(end - start + 1 for start, end in cluster)
        writer.writerow([str(i+1), str(nres), chopping])

    return output.getvalue().rstrip('\n')


def format_as_rich_table(intervals: dict, **kwargs) -> str:
    """
    Format the cluster_interval results as a rich table.
    
    Parameters:
    - intervals (dict): The cluster_interval results.
    - **kwargs: Additional keyword arguments to pass to rich.table.Table.

    Returns:
    - str: The rich formatted output.
    """
    table = Table(**kwargs)
    table.add_column("Domain", justify="right")
    table.add_column("Number of Residues", justify="right")
    table.add_column("Chopping", justify="right")

    for i, cluster in enumerate(intervals.values()):
        chopping = '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])
        nres = sum(end - start + 1 for start, end in cluster)
        table.add_row(str(i+1), str(nres), chopping)
    
    console = Console()
    with console.capture() as capture:
        console.print(table)
        console.quiet = True # This is a weird workaround to avoid rendering an empty output cell in jupyter notebooks
    console.quiet = False
    return capture.get().rstrip('\n')


def is_notebook() -> bool:
    """
    Check if the code is running in a jupyter notebook.
    
    Returns:
    - bool: True if the code is running in a jupyter notebook, False otherwise.
    """
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            return True
    except ImportError:
        return False
    return False


def format_result(intervals: dict, format: str = 'auto', delimiter: str = ',') -> str:
    """
    Format the cluster_interval results.
    Output format can be 'csv' or 'rich'. If 'auto' is selected, the output format will be 'rich' if 
    the output is a terminal or a jupyter notebook, 'csv' otherwise.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - format (str): The format of the output. Default is 'auto'. Options are 'csv', 'rich', and 'auto'. 
                    If 'auto' is selected, the output format will be 'rich' if the output is a terminal
                    "or a jupyter notebook, 'csv' otherwise.
    - delimiter (str): The delimiter to use when format is 'csv'. Default is ','.

    Returns:
    - str: The formatted output.

    Raises:
    - ValueError: If format is not one of ['auto', 'csv', 'rich'].
    """
    format_options = ['auto', 'csv', 'rich']
    format = format.lower()
    if format not in format_options:
        raise ValueError(f"format must be one of {format_options}")

    if format == 'auto':
        if is_notebook() or sys.stdout.isatty():
            return format_as_rich_table(intervals)
        else:
            return format_csv(intervals, delimiter)
        
    if format == 'csv':
        return format_csv(intervals, delimiter)
    
    if format == 'rich':
        return format_as_rich_table(intervals)
    

def print_result(intervals: dict, 
                 format: str = 'auto', 
                 delimiter: str = ',') -> None:
    """
    Print the cluster_interval results to the console.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - format (str): The format of the output. Default is 'auto'. Options are 'csv', 'rich', and 'auto'. 
                    If 'auto' is selected, the output format will be 'rich' if the output is a terminal
                    "or a jupyter notebook, 'csv' otherwise.
    - delimiter (str): The delimiter to use when format is 'csv'. Default is ','.

    Returns:
    - None
    """
    result = format_result(intervals=intervals, format=format, delimiter=delimiter)
    print(result)


def save_result(intervals: dict, 
                output_file: FilePath, 
                format: str = 'csv', 
                delimiter: str = ',') -> None:
    """
    Save the cluster_interval results to a file.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - output_file (str or Path): The file path to save the results.
    - format (str): The format of the output. Default is 'csv'. Options are 'csv', 'rich', and 'auto'. 
                    If 'auto' is selected, the output format will be 'rich' if the output is a terminal
                    "or a jupyter notebook, 'csv' otherwise.
    - delimiter (str): The delimiter to use when format is 'csv'. Default is ','.

    Returns:
    - None
    """
    with open(output_file, 'w') as f:
        # Newline at the end, else the last result line will not be shown in the file
        result = format_result(intervals=intervals, format=format, delimiter=delimiter) + '\n'
        f.write(result)
