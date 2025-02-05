from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import axes, image
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as patches


def plot_matrix(matrix: np.ndarray, 
                ax: axes.Axes = None,
                colorbar: bool = True, 
                colorbar_label: str = "", 
                colorbar_decimals: int = None, 
                **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
    """
    Plot a matrix as a heatmap.

    Parameters:
    - matrix (np.ndarray): The matrix to be plotted.
    - ax (axes.Axes, optional): The matplotlib axes object to plot the matrix on.
    - tick_top (bool, optional): Whether to place the ticks on the top of the plot.
    - colorbar (bool, optional): Whether to display a colorbar.
    - colorbar_label (str, optional): The label for the colorbar.
    - colorbar_decimals (int, optional): The number of decimal places to display in the colorbar.
    - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

    Returns:
    - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
    """
    paegreen = ["#1E461E", "#249224", "#37B137", "#56C956", "#82DD82", "#BAEFBA", "#FFFFFF"]
    cm = LinearSegmentedColormap.from_list("custom", paegreen, N=256)

    if np.all((matrix >= 0) & (matrix <= 1)):
        kwargs.setdefault("vmax", 1)
        if colorbar_decimals is None:
            colorbar_decimals = 1
    else:
        if colorbar_decimals is None:
            colorbar_decimals = 0       
    
    kwargs.setdefault("aspect", "equal")
    kwargs.setdefault("cmap", cm)
    kwargs.setdefault("vmin", 0)
    
    if ax is None:
        _, ax = plt.subplots()
    image = ax.imshow(matrix, **kwargs)

    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(image, cax=cax)
        cbar.set_label(colorbar_label, rotation=270, labelpad=15)
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter(f"%.{colorbar_decimals}f"))
    return image, ax


def plot_cluster_intervals(background_matrix: np.ndarray, 
                           cluster_intervals: dict, 
                           ax: axes.Axes = None, 
                           linewidth: int = 2, 
                           edgecolor: str = 'r', 
                           **kwargs) -> Tuple[image.AxesImage, axes.Axes]:
    """
    Plot the cluster intervals as rectangles on top of a background matrix.

    Parameters:
    - background_matrix (np.ndarray): The background matrix to plot.
    - cluster_intervals (dict): A dictionary where the keys are the cluster indices and the values are lists of
                                tuples representing the cluster intervals.
    - ax (axes.Axes, optional): The matplotlib axes object to plot the matrix on.
    - linewidth (int, optional): The width of the rectangle edges.
    - edgecolor (str, optional): The color of the rectangle edges.
    - **kwargs: Additional keyword arguments to be passed to the matplotlib.pyplot.imshow function.

    Returns:
    - Tuple[image.AxesImage, axes.Axes]: The image and axes objects.
    """
    image, ax = plot_matrix(background_matrix, ax=ax, **kwargs)
    
    for intervals in cluster_intervals.values():
        # Plot diagonal square for each interval
        for start, stop in intervals:
            rect = patches.Rectangle((start, start), stop - start, stop - start, 
                                        linewidth=linewidth, edgecolor=edgecolor, facecolor='none')
            ax.add_patch(rect)
        # Plot off-diagonal rectangles for discontinuous intervals
        for i, (start1, stop1) in enumerate(intervals[:-1]):
            for start2, stop2 in intervals[i + 1:]:
                # Off-diagonal intersection rectangles
                width1, height1 = stop1 - start1, stop2 - start2
                width2, height2 = stop2 - start2, stop1 - start1
                
                # Add the first off-diagonal rectangle
                ax.add_patch(
                    patches.Rectangle((start1, start2), width1, height1,
                                    linewidth=linewidth, edgecolor=edgecolor, facecolor='none')
                )
                # Add the second off-diagonal rectangle (symmetric to the first)
                ax.add_patch(
                    patches.Rectangle((start2, start1), width2, height2,
                                    linewidth=linewidth, edgecolor=edgecolor, facecolor='none')
                )
    return image, ax


