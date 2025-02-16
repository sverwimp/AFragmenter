
from io import StringIO
import csv
from typing import Union, Optional
from pathlib import Path
import sys

from rich.table import Table
from rich.console import Console


FilePath = Union[str, Path]


def _format_name(i: int, base_name: Optional[str] = None) -> str:
    """Formatting of the name of the cluster."""
    return '_'.join([base_name, str(i + 1)]) if base_name else str(i + 1)


def _format_nres(cluster: dict) -> int:
    """Formatting of the number of residues in the cluster."""
    return sum(end - start + 1 for start, end in cluster)


def _format_chopping(cluster: dict) -> str:
    """Formatting of the chopping of the cluster."""
    return '_'.join([f'{start + 1}-{end + 1}' for start, end in cluster])


def _row_content(cluster: dict, i: int, base_name: Optional[str] = None) -> tuple:
    """
    Format the content of a row in the result table.

    Parameters:
    - cluster (dict): The cluster information.
    - i (int): The index of the cluster.
    - base_name (str): The base name used to format the first column of the table. Default is None.

    Returns:
    - tuple: A tuple containing the name, number of residues, and chopping of the cluster.
    """
    name = _format_name(i, base_name)
    chopping = _format_chopping(cluster)
    nres = _format_nres(cluster)
    return name, str(nres), chopping


def format_csv(intervals: dict, delimiter: str = ',', base_name: Optional[str] = None) -> str:
    """
    Format the cluster_interval results as a CSV string.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - delimiter (str): The delimiter to use. Default is ','.
    - base_name (str): The base name used to format the first column of the CSV. Default is None.

    Returns:
    - str: The CSV formatted output.
    """
    output = StringIO()
    writer = csv.writer(output, delimiter=delimiter)
    writer.writerow(["domain", "nres", "chopping"])

    for i, cluster in enumerate(intervals.values()):
        name, nres, chopping = _row_content(cluster, i, base_name)
        writer.writerow([name, str(nres), chopping])

    return output.getvalue().rstrip('\n')


def _build_rich_table(intervals: dict, base_name: Optional[str] = None, **kwargs) -> str:
    """
    Format the cluster_interval results as a rich table.
    
    Parameters:
    - intervals (dict): The cluster_interval results.
    - **kwargs: Additional keyword arguments to pass to rich.table.Table.

    Returns:
    - Table: The rich Table object.
    """
    table = Table(**kwargs)
    table.add_column("Domain", justify="left")
    table.add_column("Number of Residues", justify="right")
    table.add_column("Chopping", justify="right")

    for i, cluster in enumerate(intervals.values()):
        name, nres, chopping = _row_content(cluster, i, base_name)
        table.add_row(name, str(nres), chopping)
    
    return table


def _render_rich_table(table: Table) -> str:
    """
    Render the rich table as a string.
    
    Parameters:
    - table (Table): The rich Table object.

    Returns:
    - str: The rich formatted output as a string.
    """
    console = Console()
    with console.capture() as capture:
        console.print(table)
        console.quiet = True # This is a weird workaround to avoid rendering an empty output cell in jupyter notebooks
    console.quiet = False
    return capture.get().rstrip('\n')


def format_rich_table(intervals: dict, base_name: Optional[str] = None, **kwargs) -> str:
    """
    Format the cluster_interval results as a rich table.
    
    Parameters:
    - intervals (dict): The cluster_interval results.
    - base_name (str): The base name used to format the first column of the rich table. Default is None.
    - **kwargs: Additional keyword arguments to pass to rich.table.Table.

    Returns:
    - str: The rich formatted output.
    """
    table = _build_rich_table(intervals, base_name, **kwargs)
    return _render_rich_table(table)


def is_notebook() -> bool:
    """
    Check if the code is running in a jupyter notebook.
    
    Returns:
    - bool: True if the code is running in a jupyter notebook, False otherwise.
    """
    try:
        from IPython import get_ipython # type: ignore
        if get_ipython() is not None:
            return True
    except ImportError:
        return False
    return False


def format_result(intervals: dict, format: str = 'auto', delimiter: str = ',', base_name: Optional[str] = None) -> str:
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
    
    if format == 'auto' and (is_notebook() or sys.stdout.isatty()):
        format = 'rich'
    else:
        format = 'csv'
        
    if format == 'csv':
        return format_csv(intervals, delimiter, base_name)
    
    if format == 'rich':
        return format_rich_table(intervals, base_name)
    

def print_result(intervals: dict, 
                 format: str = 'auto', 
                 delimiter: str = ',',
                 base_name: Optional[str] = None) -> None:
    """
    Print the cluster_interval results to the console.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - format (str): The format of the output. Default is 'auto'. Options are 'csv', 'rich', and 'auto'. 
                    If 'auto' is selected, the output format will be 'rich' if the output is a terminal
                    "or a jupyter notebook, 'csv' otherwise.
    - delimiter (str): The delimiter to use when format is 'csv'. Default is ','.
    - base_name (str): The base name used to format the first column of the rich table. Default is None.

    Returns:
    - None
    """
    result = format_result(intervals=intervals, format=format, delimiter=delimiter, base_name=base_name)
    print(result)


def save_result(intervals: dict, 
                output_file: FilePath, 
                base_name: Optional[str] = None,
                format: str = 'csv', 
                delimiter: str = ','
                ) -> None:
    """
    Save the cluster_interval results to a file.

    Parameters:
    - intervals (dict): The cluster_interval results.
    - output_file (str or Path): The file path to save the results.
    - base_name (str): The base name used to format the first column of the rich table. Default is None.
    - format (str): The format of the output. Default is 'csv'. Options are 'csv', 'rich', and 'auto'. 
                    If 'auto' is selected, the output format will be 'rich' if the output is a terminal
                    "or a jupyter notebook, 'csv' otherwise.
    - delimiter (str): The delimiter to use when format is 'csv'. Default is ','.

    Returns:
    - None
    """
    with open(output_file, 'w') as f:
        # Newline at the end, else the last result line will not be shown in the file
        result = format_result(intervals=intervals, format=format, delimiter=delimiter, base_name=base_name) + '\n'
        f.write(result)
