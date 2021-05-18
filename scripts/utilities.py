#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities for scripts.

Author: Andrew Tarzia

Date Created: 02 Mar 2021

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import re
from scipy.spatial.distance import euclidean


def colors_i_like(palette=None):
    """
    A list of colours I like to choose from.

    palette options:
        None
        Base
        IBM
        Wong
        Tol
        CB_pairs

    """

    if palette is None:
        return [
            '#FA7268', '#F8A72A', '#DAF7A6', '#900C3F', '#6BADB0',
            '#DB869D', '#F6D973', 'mediumvioletred',
            'skyblue', 'gold', 'palegreen', 'coral',
        ]
    elif palette == 'Base':
        return [
            '#D81B60', '#1E88E5', '#FFC107', '#FE6100', '#004D40'
        ]
    elif palette == 'IBM':
        return [
            '#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000'
        ]
    elif palette == 'Wong':
        return [
            '#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442',
            '#0072B2', '#D55E00', '#CC79A7'
        ]
    elif palette == 'Tol':
        return [
            '#332288', '#117733', '#44AA99', '#88CCEE', '#DDCC77',
            '#CC6677', '#AA4499', '#882255',
        ]
    elif palette == 'CB_pairs':
        return [
            '#FFC20A', '#0C7BDC', '#994F00', '#006CD1', '#E1BE6A',
            '#40B0A6', '#E66100', '#5D3A9B', '#1AFF1A', '#4B0092',
            '#FEFE62', '#D35FB7', '#005AB5', '#DC3220', '#1A85FF',
            '#D41159',
        ]


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249
    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, normal=None):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    If normal is given, the angle polarity is determined using the
    cross product of the two vectors.

    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if normal is not None:
        # Get normal vector and cross product to determine sign.
        cross = np.cross(v1_u, v2_u)
        if np.dot(normal, cross) < 0:
            angle = -angle
    return angle


def get_atom_distance(molecule, atom1_id, atom2_id):
    """
    Return the distance between atom1 and atom2.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`

    atom1_id : :class:`int`
        The id of atom1.

    atom2_id : :class:`int`
        The id of atom2.

    Returns
    -------
    :class:`float`
        The euclidean distance between two atoms.

    """

    position_matrix = molecule.get_position_matrix()

    distance = euclidean(
        u=position_matrix[atom1_id],
        v=position_matrix[atom2_id]
    )

    return float(distance)


def replace(string, substitutions):
    """
    Replace multiple substrings.

    Parameters
    ----------
    string : :class:`str`
        String to make substitutions in.

    substitutions : :class:`dict` of :class:`str` pairs
        Dictionary where values replace keys in string.

    Returns
    -------
    :class:`str`
        String with substitutions made.

    """

    substrings = sorted(substitutions, key=len, reverse=True)
    regex = re.compile('|'.join(map(re.escape, substrings)))
    return regex.sub(
        lambda match: substitutions[match.group(0)], string
    )


def bar_chart(energies, filename, x_labels, y_title, facecolor):

    xs = [1, 2, 3, 4]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(
        xs,
        energies,
        width=1,
        color=facecolor,
        edgecolor='k',
        linewidth=2
    )
    ax.set_xticks(xs)
    ax.set_xticklabels(x_labels)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand', fontsize=16)
    ax.set_ylabel(y_title, fontsize=16)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
    plt.close()


def get_g16_energy(
    filename,
    front_splitter,
    back_splitter,
    back_splitter2,
):

    with open(filename, 'r') as f:
        data = f.read()

    energy = data.split(front_splitter)
    energy = energy[-1].split(back_splitter)[0]
    try:
        return float(energy)  # a.u.
    except ValueError:
        if back_splitter2 is None:
            raise ValueError()
        energy = energy.split(back_splitter2)[0]
        return float(energy)  # a.u.


def collate_energies(file_list, functional):

    collated_energies = {}
    for file in file_list:
        splits = file.split('/')
        mol_name = splits[1]
        splits = mol_name.split('_')
        mol_name = splits[1].lower()+'_'+splits[2][0].upper()
        energy = get_g16_energy(
            filename=file,
            front_splitter=f'SCF Done:  E({functional}) =',
            back_splitter='A.U. after',
            back_splitter2=None,
        )
        collated_energies[mol_name] = energy

    return collated_energies


def scatter_plot(
    X,
    Y,
    xtitle,
    ytitle,
    xlim,
    ylim,
    title=None,
    c='firebrick',
    edgecolors='k',
    marker='o',
    alpha=1.0,
    s=80,
    Z=None,
    cmap=None
):
    """
    Make scatter plot.

    Parameters
    ----------
    X : :class:``

    Y : :class:``

    xtitle : :class:``

    ytitle : :class:``

    xlim : :class:``

    ylim : :class:``

    title : :class:``

    c : :class:``

    edgecolors : :class:``

    marker : :class:``

    alpha : :class:``

    s : :class:``

    Z : :class:``

    cmap : :class:`dict`
        Dictionary containing information for cmap.
        Example:
        {
            'mid_point': 0.5,
            'cmap': cm.Purples,
            'ticks': [0, .50, 1.00],
            'labels': [
                '0',
                '20',
                '40'
            ],
            'cmap_label': 'flex',
        }

        This requries that the `Z` argument is on the range of 0 to 1.
        Use 'labels' to define the relationship between cmap and
        Z value.


    Returns
    -------

    fig

    ax

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    if cmap is None and Z is None:
        ax.scatter(
            X, Y,
            c=c,
            edgecolors=edgecolors,
            marker=marker,
            alpha=alpha,
            s=s
        )
    else:
        cmp = define_plot_cmap(
            fig, ax,
            mid_point=cmap['mid_point'],
            cmap=cmap['cmap'],
            ticks=cmap['ticks'],
            labels=cmap['labels'],
            cmap_label=cmap['cmap_label']
        )
        ax.scatter(
            X, Y,
            c=cmp(Z),
            edgecolors=edgecolors,
            marker=marker,
            alpha=alpha,
            s=s
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if title is not None:
        ax.set_title(title, fontsize=16)
    return fig, ax


def define_plot_cmap(fig, ax, mid_point, cmap, ticks, labels,
                     cmap_label):
    """
    Define cmap shifted to midpoint and plot colourbar

    """
    new_cmap = shiftedColorMap(
        cmap,
        midpoint=mid_point,
        name='shifted'
    )
    X = np.linspace(0, 1, 256)
    cax = ax.scatter(-X-100, -X-100, c=X, cmap=new_cmap)
    cbar = fig.colorbar(cax, ticks=ticks, spacing='proportional')
    cbar.ax.set_yticklabels(labels, fontsize=16)
    cbar.set_label(cmap_label, fontsize=16)
    return new_cmap


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0,
                    name='shiftedcmap'):
    """
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    From Stack Exchange:
        https://stackoverflow.com/questions/7404116/
        defining-the-midpoint-of-a-colormap-in-matplotlib

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    """
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)
    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)
    return newcmap
