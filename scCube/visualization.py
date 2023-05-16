import operator

import matplotlib.colors as clr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pandas import DataFrame
from typing import Optional, Tuple, Sequence, Union
import seaborn as sns


def plot_spatial_pattern_scatter(
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        dim: int = 2,
        x: str = "point_x",
        y: str = "point_y",
        z: Optional[str] = None,
        label: Optional[str] = None,
        palette: Optional[list] = None,
        colormap: str = 'rainbow',
        size: float = 10,
        alpha: float = 1,
):
    """
    Plot scatter plot of cell type spatial pattern
    :param obj: DataFrame of meta
    :param figwidth: Figure width
    :param figheight: Figure height
    :param dim: Spatial dimensionality
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param z: The name of column containing z coordinate, only use when 'dim = 3'
    :param label: The name of column containing cell type information, if 'label == None', plot coordinates without
    cell type information only
    :param palette: List of colors used, if 'palette == None', plot scatter plot with colormap colors
    :param colormap: The name of cmap
    :param size: The size of point
    :param alpha: The transparency of point
    :return: fig
    """

    if label is None:
        if dim == 2:
            fig, ax = plt.subplots(figsize=(figwidth, figheight))
            ax.scatter(
                obj[x],
                obj[y],
                s=size,
                alpha=alpha)
            ax.set_xlabel(x)
            ax.set_ylabel(y)
        elif dim == 3:
            fig = plt.figure(figsize=(figwidth, figheight))
            ax = fig.add_subplot(projection='3d')
            ax.scatter(
                obj[x],
                obj[y],
                obj[z],
                s=size,
                alpha=alpha)
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            ax.set_zlabel(z)
        plt.title("Spatial coordinates of generated data")

    else:
        label_list = sorted(list(obj[label].unique()))
        if palette is None:
            color_map = plt.cm.get_cmap(colormap, len(label_list))
            palette = [clr.rgb2hex(color_map(i)) for i in range(color_map.N)]
        else:
            assert len(palette) >= len(label_list), "The number of colors provided is too few!!"

        if dim == 2:
            fig, ax = plt.subplots(figsize=(figwidth, figheight))
            for i in range(len(label_list)):
                ax.scatter(
                    obj.loc[obj[label] == label_list[i], x],
                    obj.loc[obj[label] == label_list[i], y],
                    color=palette[i],
                    label=label_list[i],
                    s=size,
                    alpha=alpha)

            ax.set_xlabel(x)
            ax.set_ylabel(y)

        elif dim == 3:
            fig = plt.figure(figsize=(figwidth, figheight))
            ax = fig.add_subplot(projection='3d')
            for i in range(len(label_list)):
                ax.scatter(
                    obj.loc[obj[label] == label_list[i], x],
                    obj.loc[obj[label] == label_list[i], y],
                    obj.loc[obj[label] == label_list[i], z],
                    color=palette[i],
                    label=label_list[i],
                    s=size,
                    alpha=alpha)

            ax.set_xlabel(x)
            ax.set_ylabel(y)
            ax.set_zlabel(z)

        plt.title("Spatial pattern of generated data")
        plt.legend(bbox_to_anchor=(1, 0.6), borderaxespad=0, frameon=False)

    return fig


def plot_spatial_pattern_density(
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        x: str = "point_x",
        y: str = "point_y",
        show_celltype: Optional[str] = None,
        label: str = "Cell_type",
        colormap: str = 'Blues',
        fill: bool = True,
):
    """
    Plot density plot of cell type spatial pattern
    :param obj: DataFrame of meta
    :param figwidth: Figure width
    :param figheight: Figure height
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param show_celltype: The cell type selected to plot separately, if 'show_celltype == None',
    plot all cell type together
    :param label: The name of column containing cell type information
    :param colormap: The name of cmap
    :param fill: If 'fill = True', fill in the area between bivariate contours
    :return: fig
    """

    label_list = sorted(list(obj[label].unique()))
    if show_celltype is None:
        fig, axes = plt.subplots(2, int(np.ceil(len(label_list) / 2)), figsize=(figwidth, figheight))
        palette = sns.color_palette(colormap, len(label_list))
        for ax, s, ct in zip(axes.flat, palette, label_list):
            # Create a cubehelix colormap to use with kdeplot
            # colormap = sns.cubehelix_palette(start=s, light=1, as_cmap=True)
            colormap = sns.light_palette(s, as_cmap=True)

            x_new = obj.loc[obj[label] == ct, x]
            y_new = obj.loc[obj[label] == ct, y]

            sns.kdeplot(
                x=x_new,
                y=y_new,
                cmap=colormap,
                fill=fill,
                ax=ax,
            )
            ax.set_title(ct)
        plt.tight_layout()
        plt.suptitle('The density of each cell type', fontsize=16, y=1)

    else:
        assert show_celltype in label_list, "The selected cell type does not exist!!"
        fig = sns.jointplot(
            data=obj.loc[obj[label] == show_celltype, ],
            x=x,
            y=y,
            kind='kde',
            cmap=colormap,
            fill=fill,
            height=figheight,
            width=figwidth)

        fig.fig.suptitle('The density of ' + show_celltype, fontsize=16, y=1.05)

    return fig


def get_palette(categories, cmap):
    """
    Generate dictionary mapping categories to color.
    """
    cc = plt.cm.get_cmap(cmap, len(categories))
    colors = [clr.rgb2hex(cc(i)) for i in range(cc.N)]
    cmap_new = clr.ListedColormap(colors)
    palette = {x: cmap_new(i) for i, x in enumerate(categories)}
    return palette


# https://github.com/racng/scatterpie
def pie_marker(
        ratios: Sequence[float],
        res: int = 50,
        direction: str = "+",
        start: float = 0.0,
) -> Tuple[list, list]:
    """
    Create each slice of pie as a separate marker.
    Parameters:
        ratios(list): List of ratios that add up to 1.
        res: Number of points around the circle.
        direction: '+' for counter-clockwise, or '-' for clockwise.
        start: Starting position in radians.
    Returns:
        xys, ss: Tuple of list of xy points and sizes of each slice in the pie marker.
    """

    if np.abs(np.sum(ratios) - 1) > 0.01:
        print("Warning: Ratios do not add up to 1.")

    if direction == '+':
        op = operator.add
    elif direction == '-':
        op = operator.sub

    xys = []  # list of xy points of each slice
    ss = []  # list of size of each slice
    start = float(start)
    for ratio in ratios:
        # points on the circle including the origin (0,0) and the slice
        end = op(start, 2 * np.pi * ratio)
        n = round(ratio * res)  # number of points forming the arc
        x = [0] + np.cos(np.linspace(start, end, n)).tolist()
        y = [0] + np.sin(np.linspace(start, end, n)).tolist()
        xy = np.column_stack([x, y])
        xys.append(xy)
        ss.append(np.abs(xy).max())
        start = end

    return xys, ss


def plot_spot_scatterpie(
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        x: str = "spot_x",
        y: str = "spot_y",
        # cols: Optional[list] = None,
        palette: Optional[dict] = None,
        colormap: str = 'rainbow',
        res: int = 50,
        direction: str = "+",
        start: float = 0.0,
        size=100,
        edgecolor="none",
):
    """
    Plot scatterpie plot of spot-based data
    :param obj: DataFrame of cell type proportion per spot
    :param figwidth: Figure width
    :param figheight: Figure height
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param palette: Dict of color of each cell type, if 'palette == None', plot scatterpie plot with colormap colors
    :param colormap: The name of cmap
    :param res: Number of points around the circle
    :param direction: '+' for counter-clockwise, or '-' for clockwise
    :param start: Starting position in radians
    :param size: The size of point
    :param edgecolor: The edge color of point
    :return: ax
    """

    # make copy of dataframe and set xy as index
    obj = obj.copy().set_index([x, y])
    label_list = obj.columns
    obj = obj.reset_index()

    if palette is None:
        palette = get_palette(label_list, colormap)

    ratios = obj[label_list].to_records(index=False).tolist()
    colors = [palette[cat] for cat in label_list]

    _, ax = plt.subplots(figsize=(figwidth, figheight))
    # make pie marker for each unique set of ratios
    df = pd.DataFrame({'x': obj[x].values, 'y': obj[y].values, 'ratios': ratios})
    df.ratios = df.ratios.apply(tuple)
    gb = df.groupby("ratios")
    for ratio in gb.groups:
        group = gb.get_group(ratio)
        xys, ss = pie_marker(ratio, res=res, direction=direction, start=start)
        for xy, s, color in zip(xys, ss, colors):
            # plot non-zero slices
            if s != 0:
                ax.scatter(group.x, group.y, marker=xy, s=[s * s * size], facecolor=color, edgecolor=edgecolor)

    handles = [plt.scatter([], [], color=palette[i], label=i) for i in label_list]
    ax.legend(handles=handles, bbox_to_anchor=(1, 0.6), borderaxespad=0, frameon=False)

    return ax


def plot_spot_prop(
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        x: str = "spot_x",
        y: str = "spot_y",
        colormap: str = 'viridis',
        show_celltype: Union[list, str] = "",
        size=100,
        alpha=1,
):
    """
    Plot scatter plot of proportion of selected cell type
    :param obj: DataFrame of cell type proportion per spot
    :param figwidth: Figure width
    :param figheight: Figure height
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param colormap: The name of cmap
    :param show_celltype: The cell type selected to plot
    :param size: The size of point
    :param alpha: The transparency of point
    :return: fig
    """

    x_new = obj[x]
    y_new = obj[y]
    if type(show_celltype) == str:
        cell_type = obj[show_celltype]
        fig, ax = plt.subplots(figsize=(figwidth, figheight))
        g = ax.scatter(x_new, y_new, s=size, cmap=colormap, c=cell_type, alpha=alpha)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(show_celltype)
        fig.colorbar(g)

    elif type(show_celltype) == list:
        cell_type_list = sorted(show_celltype)

        fig, axes = plt.subplots(1, len(cell_type_list), figsize=(figwidth, figheight))

        for i in range(len(cell_type_list)):
            cell_type_tmp = obj[cell_type_list[i]]
            g = axes[i].scatter(x_new, y_new, s=size, cmap=colormap, c=cell_type_tmp, alpha=alpha)
            axes[i].set_xlabel(x)
            axes[i].set_ylabel(y)
            axes[i].set_title(cell_type_list[i])
            fig.colorbar(g, ax=axes[i])

    plt.tight_layout()

    return g


def plot_gene_scatter(
        data: DataFrame,
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        dim: int = 2,
        label: str = 'Cell',
        normalize: bool = True,
        x: str = "point_x",
        y: str = "point_y",
        z: str = "point_z",
        show_gene: str = "",
        colormap: str = 'viridis',
        size: float = 10,
        alpha: float = 1
):
    """
    Plot scatter plot of spatial expression pattern of selected gene
    :param data: DataFrame of data
    :param obj: DataFrame of meta
    :param figwidth: Figure width
    :param figheight: Figure height
    :param dim: Spatial dimensionality
    :param label: The name of column containing cell/spot name
    :param normalize: If 'normalize = True', normalizing expression value to [0, 1]
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param z: The name of column containing z coordinate
    :param show_gene: The gene selected to plot
    :param colormap: The name of cmap
    :param size: The size of point
    :param alpha: The transparency of point
    :return: fig
    """

    obj.index = list(obj[label])
    gene_exp = np.array(data.T[show_gene])

    if normalize:
        gene_exp = (gene_exp - gene_exp.min()) / (gene_exp.max() - gene_exp.min())

    if dim == 2:
        fig, ax = plt.subplots(figsize=(figwidth, figheight))

        x_new = obj[x]
        y_new = obj[y]
        g = ax.scatter(x_new, y_new, s=size, cmap=colormap, c=gene_exp, alpha=alpha)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        fig.colorbar(g)
    elif dim == 3:
        fig = plt.figure(figsize=(figwidth, figheight))
        ax = fig.add_subplot(projection='3d')

        x_new = obj[x]
        y_new = obj[y]
        z_new = obj[z]
        g = ax.scatter(x_new, y_new, z_new, s=size, cmap=colormap, c=gene_exp, alpha=alpha)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_zlabel(z)
        pos = ax.get_position()
        ax2 = fig.add_axes([pos.xmax + 0.25, pos.ymin, 0.03, 0.7 * (pos.ymax - pos.ymin)])
        fig.colorbar(g, cax=ax2)

    plt.title(show_gene)
    plt.tight_layout()

    return g


def plot_spot_histplot(
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        label: str = 'spot',
        n_bins: int = 20
):
    """
    Plot histplot of spot-based data to investigate cell number per spot
    :param obj: DataFrame of cell-spot index
    :param figwidth: Figure width
    :param figheight: Figure height
    :param label: The name of column containing spot name
    :param n_bins: The number of equal-width bins in the range
    :return: fig
    """

    x = np.array(obj[label].value_counts())
    mu = np.mean(x)
    sigma = np.std(x, ddof=1)

    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    # the histogram of the data
    n, bins, patches = ax.hist(x, n_bins, density=True)
    # add a 'best fit' line
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu)) ** 2))
    ax.plot(bins, y, '--')
    ax.set_xlabel('Cell number per spot')
    ax.set_ylabel('Density')
    ax.set_title(r'$\mu='+ str(round(mu, 2)) + '$, $\sigma='+ str(round(sigma, 2)) + '$')

    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()

    return fig


def plot_slice_scatter(
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        x: str = "point_x",
        y: str = "point_y",
        label: str = 'Cell_type',
        palette: Optional[list] = None,
        colormap: str = 'rainbow',
        size: float = 10,
        alpha: float = 1
):
    """
    Plot 2d scatter plot of cell type spatial pattern of each slices from 3d data
    :param obj: DataFrame of meta
    :param figwidth: Figure width
    :param figheight: Figure height
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param label: The name of column containing cell type information
    :param palette: List of colors used, if 'palette == None', plot scatter plot with colormap colors
    :param colormap: The name of cmap
    :param size: The size of point
    :param alpha: The transparency of point
    :return: fig
    """
    slice_id = sorted(list(np.unique(obj['slice'].values)))
    label_list = sorted(list(np.unique(obj[label].values)))
    if palette is None:
        color_map = plt.cm.get_cmap(colormap, len(label_list))
        palette = [clr.rgb2hex(color_map(i)) for i in range(color_map.N)]
    else:
        assert len(palette) >= len(label_list), "The number of colors provided is too few!!"

    fig, axes = plt.subplots(1, len(slice_id), figsize=(figwidth, figheight))

    for i in range(len(slice_id)):
        for j in range(len(label_list)):
            axes[i].scatter(
                obj.loc[(obj['slice'] == slice_id[i]) & (obj[label] == label_list[j]), x],
                obj.loc[(obj['slice'] == slice_id[i]) & (obj[label] == label_list[j]), y],
                color=palette[j],
                label=label_list[j],
                s=size,
                alpha=alpha)

            axes[i].set_xlabel(x)
            axes[i].set_ylabel(y)
            axes[i].set_title('Slice ' + str(slice_id[i]))

    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1, 0.6), borderaxespad=0, frameon=False)

    return fig


def plot_slice_gene_scatter(
        data: DataFrame,
        obj: DataFrame,
        figwidth: float = 8,
        figheight: float = 8,
        x: str = "point_x",
        y: str = "point_y",
        label: str = 'Cell_type',
        normalize: bool = True,
        show_gene: str = "",
        colormap: str = 'viridis',
        size: float = 10,
        alpha: float = 1
):
    """
    Plot 2d scatter plot of spatial expression pattern of selected gene of each slices from 3d data
    :param data: DataFrame of data
    :param obj: DataFrame of meta
    :param figwidth: Figure width
    :param figheight: Figure height
    :param x: The name of column containing x coordinate
    :param y: The name of column containing y coordinate
    :param label: The name of column containing cell type information
    :param normalize: If 'normalize = True', normalizing expression value to [0, 1]
    :param show_gene: The gene selected to plot
    :param colormap: The name of cmap
    :param size: The size of point
    :param alpha: The transparency of point
    :return: fig
    """
    slice_id = sorted(list(np.unique(obj['slice'].values)))
    fig, axes = plt.subplots(1, len(slice_id), figsize=(figwidth, figheight))

    for i in range(len(slice_id)):
        obj_new = obj.loc[obj['slice'] == slice_id[i], ]
        data_new = data.T.iloc[obj_new.index].T
        gene_exp = np.array(data_new.T[show_gene])

        if normalize:
            gene_exp = (gene_exp - gene_exp.min()) / (gene_exp.max() - gene_exp.min())

        x_new = obj_new[x]
        y_new = obj_new[y]
        g = axes[i].scatter(x_new, y_new, s=size, cmap=colormap, c=gene_exp, alpha=alpha)
        axes[i].set_xlabel(x)
        axes[i].set_ylabel(y)
        axes[i].set_title(show_gene + ' (Slice ' + str(slice_id[i]) + ')')
        fig.colorbar(g, ax=axes[i])

    plt.tight_layout()

    return g


