#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting cage properties.

Author: Andrew Tarzia

Date Created: 29 May 2020

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from pandas import read_csv


def define_plot_cmap(
    fig,
    ax,
    mid_point,
    cmap,
    ticks,
    labels,
    cmap_label,
):
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


def shiftedColorMap(
    cmap,
    start=0,
    midpoint=0.5,
    stop=1.0,
    name='shiftedcmap',
):
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


def histogram_plot_N(
    Y,
    X_range,
    width,
    alpha,
    color,
    edgecolor,
    xtitle,
    labels=None,
    density=False,
    N=1
):
    """
    Make histogram plot with 1 distribution.

    """

    fig, ax = plt.subplots(figsize=(8, 5))
    X_bins = np.arange(X_range[0], X_range[1], width)
    if N == 1:
        hist, bin_edges = np.histogram(
            a=Y,
            bins=X_bins,
            density=density
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            align='edge',
            alpha=alpha,
            width=width,
            color=color,
            edgecolor=edgecolor
        )
    else:
        for i_ in range(N):
            if type(color) is not list or len(Y) != N:
                raise ValueError(
                    'Make sure color and Y are of length N'
                )
            hist, bin_edges = np.histogram(
                a=Y[i_],
                bins=X_bins,
                density=density
            )
            if labels[i_] is None:
                label = ''
            else:
                label = labels[i_]
            ax.bar(
                bin_edges[:-1],
                hist,
                align='edge',
                alpha=alpha[i_],
                width=width,
                color=color[i_],
                edgecolor=edgecolor,
                label=label
            )
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    if density is False:
        ax.set_ylabel('count', fontsize=16)
    elif density is True:
        ax.set_ylabel('frequency', fontsize=16)
    ax.set_xlim(X_range)
    if N > 1 and labels[0] is not None:
        ax.legend(fontsize=16)
    return fig, ax


def plot_energetics_and_geom(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    energy_preferences,
    plane_devs,
    sqpl_ops,
):
    """
    Plot energy preference of all cages.

    """

    names = {
        'plane_dev': {
            'xlim': (0, 1),
            'xtitle': r'$D_{\mathrm{max}}$ [$\mathrm{\AA}$]',
        },
        'sqpl': {
            'xlim': (0.6, 1),
            'xtitle': r'$q_{\mathrm{sqp,min}}$',
        }
    }

    for name, data in zip(names, [plane_devs, sqpl_ops]):
        fig, ax = plt.subplots(figsize=(8, 5))

        c_passed = colors_i_like()[4]
        # c_failed = colors_i_like()[3]
        c_negative = colors_i_like()[3]
        c_experiments = colors_i_like()[0]
        m_passed = 'o'
        # m_failed = 'o'
        m_negative = 's'
        m_experiments = 'X'

        x_passed = []
        y_passed = []
        # x_failed = []
        # y_failed = []
        x_negative = []
        y_negative = []
        x_experiments = []
        y_experiments = []

        for i, lig in enumerate(ligands):
            if lig in experiments:
                x_experiments.append(data[i])
                y_experiments.append(energy_preferences[i])
            elif energy_preferences[i] < 0:
                x_negative.append(data[i])
                y_negative.append(energy_preferences[i])
            elif lig in cages_cis_wins or lig in cages_not_wins:
                x_passed.append(data[i])
                y_passed.append(energy_preferences[i])
            # elif lig in cages_not_wins:
            #     x_failed.append(data[i])
            #     y_failed.append(energy_preferences[i])
            else:
                raise ValueError('no matches!?')

        print(
            sum([len(x_experiments), len(x_negative), len(x_passed)])
        )

        ax.scatter(
            x_passed,
            y_passed,
            c=c_passed,
            edgecolors='k',
            marker=m_passed,
            alpha=1,
            s=80,
            label='$cis$ preferred'
        )
        ax.scatter(
            x_negative,
            y_negative,
            c=c_negative,
            edgecolors='k',
            marker=m_negative,
            alpha=1,
            s=80,
            label='$cis$ not preferred'
        )
        ax.scatter(
            x_experiments,
            y_experiments,
            c=c_experiments,
            edgecolors='k',
            marker=m_experiments,
            alpha=1,
            s=80,
            label='published examples'
        )

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(names[name]['xtitle'], fontsize=16)
        ax.set_ylabel('stability of C isomer [kJmol$^{-1}$]', fontsize=16)
        # ax.set_xlim(names[name]['xlim'])
        # ax.set_ylim(-40, 80)

        ax.axhline(y=6.0, c='k', alpha=0.6, lw=2)
        # if name == 'sqpl':
        #     ax.axvline(x=0.95, c='k', alpha=0.6, lw=2)
        ax.legend(fontsize=16)

        fig.tight_layout()
        fig.savefig(
            f'all_cages_pref_and_stable_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def plot_energetics_and_geom_3D(
    ligands,
    experiments,
    energy_preferences,
    plane_devs,
    sqpl_ops,
):
    """
    Plot energy preference of all cages.

    """

    names = {
        'plane_dev': {
            'xlim': (0, None),
            'xtitle': r'$D_{\mathrm{max}}$ [$\mathrm{\AA}$]',
        },
        'sqpl': {
            'xlim': (None, 1),
            'xtitle': r'$q_{\mathrm{sqp,min}}$',
        }
    }

    m_1 = 'o'
    m_2 = 'X'

    x_1 = []
    y_1 = []
    z_1 = []
    x_2 = []
    y_2 = []
    z_2 = []

    for i, lig in enumerate(ligands):
        if lig in experiments:
            x_2.append(sqpl_ops[i])
            y_2.append(energy_preferences[i])
            z_2.append(plane_devs[i])
        else:
            x_1.append(sqpl_ops[i])
            y_1.append(energy_preferences[i])
            z_1.append(plane_devs[i])

    print(sum([len(x_1), len(x_2)]))

    # Normalize plane deviations.
    z_1 = [i/0.1 for i in z_1]
    z_2 = [i/0.1 for i in z_2]

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = {
        'mid_point': 0.5,
        'cmap': cm.Blues_r,
        'ticks': [0, .50, 1.00],
        'labels': ['0', '0.05', '0.1'],
        'cmap_label': names['plane_dev']['xtitle'],
    }
    cmp = define_plot_cmap(
        fig, ax,
        mid_point=cmap['mid_point'],
        cmap=cmap['cmap'],
        ticks=cmap['ticks'],
        labels=cmap['labels'],
        cmap_label=cmap['cmap_label']
    )

    ax.scatter(
        x_1,
        y_1,
        c=cmp(z_1),
        edgecolors='k',
        marker=m_1,
        alpha=1,
        s=140,
        # label='$cis$ preferred'
    )
    ax.scatter(
        x_2,
        y_2,
        c=cmp(z_2),
        edgecolors='k',
        marker=m_2,
        alpha=1,
        s=140,
        # label='published examples'
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(names['sqpl']['xtitle'], fontsize=16)
    # ax.set_xlabel(names[name]['xtitle'], fontsize=16)
    ax.set_ylabel('stability of C isomer [kJmol$^{-1}$]', fontsize=16)
    ax.set_xlim(0.3, 1.05)
    ax.set_ylim(-40, 50)

    ax.axhline(y=6.0, c='r', alpha=0.8, lw=2)
    ax.axhline(y=0.0, c='k', alpha=0.8, lw=2, linestyle='--')

    ax.scatter(
        -100, -100,
        c='white',
        edgecolors='k',
        marker=m_1,
        alpha=1,
        s=140,
        label='this work',
    )
    ax.scatter(
        -100, -100,
        c='white',
        edgecolors='k',
        marker=m_2,
        alpha=1,
        s=140,
        label='published examples',
    )
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        'all_cages_pref_and_stable_3D.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = {
        'mid_point': 0.5,
        'cmap': cm.Blues_r,
        'ticks': [0, .50, 1.00],
        'labels': ['0', '0.05', '0.1'],
        'cmap_label': names['plane_dev']['xtitle'],
    }
    cmp = define_plot_cmap(
        fig, ax,
        mid_point=cmap['mid_point'],
        cmap=cmap['cmap'],
        ticks=cmap['ticks'],
        labels=cmap['labels'],
        cmap_label=cmap['cmap_label']
    )

    ax.scatter(
        x_1,
        y_1,
        c=cmp(z_1),
        edgecolors='k',
        marker=m_1,
        alpha=1,
        s=140,
        # label='$cis$ preferred'
    )
    ax.scatter(
        x_2,
        y_2,
        c=cmp(z_2),
        edgecolors='k',
        marker=m_2,
        alpha=1,
        s=140,
        # label='published examples'
    )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(names['sqpl']['xtitle'], fontsize=16)
    # ax.set_xlabel(names[name]['xtitle'], fontsize=16)
    # ax.set_ylabel(
    #    'stability of C isomer [kJmol$^{-1}$]',
    #    fontsize=16
    # )
    ax.set_ylabel(
        r'$\Delta E_{\mathrm{cis}}$ [kJmol$^{-1}$]',
        fontsize=16,
    )
    ax.set_xlim(0.95, 1)
    ax.set_ylim(-20, 30)

    # ax.axhline(y=6.0, c='r', alpha=0.8, lw=2)
    # ax.axhline(y=0.0, c='k', alpha=0.8, lw=2, linestyle='--')

    fig.tight_layout()
    fig.savefig(
        'all_cages_pref_and_stable_3D_zoomed.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_all_cages_bars(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    y_value,
    y_title,
    y_bar,
    suffix,
):
    """
    Plot y value as bar chart of all cages.

    """

    fig, ax = plt.subplots(figsize=(8, 5))

    c_passed = colors_i_like()[4]
    c_failed = colors_i_like()[3]
    c_experiments = colors_i_like()[0]

    x_passed = []
    y_passed = []
    x_failed = []
    y_failed = []
    x_experiments = []
    y_experiments = []

    for i, lig in enumerate(ligands):
        if lig in experiments:
            x_experiments.append(i+1)
            y_experiments.append(y_value[i])
        elif lig in cages_cis_wins:
            x_passed.append(i+1)
            y_passed.append(y_value[i])
        elif lig in cages_not_wins:
            x_failed.append(i+1)
            y_failed.append(y_value[i])
        else:
            raise ValueError('no matches!?')

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cage', fontsize=16)
    ax.set_ylabel(y_title, fontsize=16)
    ax.set_xlim(0, 61)

    ax.bar(
        x_passed,
        y_passed,
        color=c_passed,
        width=1,
        edgecolor=c_passed,
        alpha=1,
        label='C isomer passed'
    )
    ax.bar(
        x_failed,
        y_failed,
        color=c_failed,
        width=1,
        edgecolor=c_failed,
        alpha=1,
        label='C isomer failed'
    )
    ax.bar(
        x_experiments,
        y_experiments,
        color=c_experiments,
        width=1,
        edgecolor=c_experiments,
        alpha=1,
        label='published examples'
    )
    ax.legend(fontsize=16)
    ax.axhline(y=y_bar, c='k', alpha=0.8, lw=2)

    fig.tight_layout()
    fig.savefig(
        f'all_C_cages_bar_{suffix}.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def isomer_plot(dictionary, file_name, ytitle, ylim, horiz=None):
    """
    Generic plot of isomer properties.

    """

    col = 3

    X_positions = {'A': 2, 'B': 4, 'C': 6, 'D': 8}

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        list(X_positions.values()),
        list(dictionary.values()),
        color=colors_i_like()[col],
        lw=3,
    )
    for isomer in dictionary:
        Y = dictionary[isomer]
        ax.scatter(
            X_positions[isomer],
            Y,
            c=colors_i_like()[col],
            edgecolors='none',
            marker='o',
            alpha=1,
            s=180
        )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xticks([X_positions[i] for i in X_positions])
    ax.set_xticklabels(list(X_positions.keys()))
    ax.set_xlim(0, 10)
    ax.set_ylim(ylim)
    if horiz is not None:
        for i, j in zip(*horiz):
            ax.axhline(y=i, c=j[0], lw=j[1], alpha=0.4)
    fig.tight_layout()
    fig.savefig(
        file_name,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_isomer_distributions():
    """
    Plot y value as bar chart of all cages.

    """

    data = read_csv('all_cage_results.txt')

    names = {
        'plane_dev': {
            'xlim': (0, 1),
            'xtitle': r'$D_{\mathrm{max}}$ [$\mathrm{\AA}$]',
            'atitle': 'plane_dev_A',
            'btitle': 'plane_dev_B',
            'ctitle': 'plane_dev_C',
            'dtitle': 'plane_dev_D',
            'width': 0.05,
        },
        'sqpl': {
            'xlim': (0.0, 1),
            'xtitle': r'$q_{\mathrm{sqp,min}}$',
            'atitle': 'sqpl_op_A',
            'btitle': 'sqpl_op_B',
            'ctitle': 'sqpl_op_C',
            'dtitle': 'sqpl_op_D',
            'width': 0.05,
        },
        'energy': {
            'xlim': (0, 200),
            'xtitle': r'relative isomer energy [kJmol$^{-1}$]',
            'atitle': 'energy_A',
            'btitle': 'energy_B',
            'ctitle': 'energy_C',
            'dtitle': 'energy_D',
            'width': 5,
        },
    }

    for name in names:
        fig, ax = plt.subplots(figsize=(5, 5))
        # X_bins = np.arange(
        #     names[name]['xlim'][0],
        #     names[name]['xlim'][1],
        #     names[name]['width'],
        # )

        c_a = colors_i_like()[4]
        c_b = colors_i_like()[3]
        c_c = colors_i_like()[0]
        c_d = colors_i_like()[2]
        x_a = list(data[names[name]['atitle']])
        x_b = list(data[names[name]['btitle']])
        x_c = list(data[names[name]['ctitle']])
        x_d = list(data[names[name]['dtitle']])
        labels = ['A', 'B', 'C', 'D']
        XS = [1, 2, 3, 4]
        YS = [x_a, x_b, x_c, x_d]
        CS = [c_a, c_b, c_c, c_d]

        for X, label, Y, C in zip(XS, labels, YS, CS):
            # hist, bin_edges = np.histogram(
            #     a=Y,
            #     bins=X_bins,
            #     density=True
            # )
            # ax.bar(
            #     bin_edges[:-1],
            #     hist,
            #     align='edge',
            #     alpha=0.2,
            #     width=names[name]['width'],
            #     color=C,
            #     edgecolor='none',
            #     # linewidth=2,
            #     label=label,
            # )
            # ax.plot(
            #     bin_edges[:-1],
            #     hist,
            #     # align='edge',
            #     alpha=1.0,
            #     # width=names[name]['width'],
            #     color=C,
            #     # edgecolor=C,
            #     marker='o',
            #     linewidth=2,
            #     label=label,
            # )
            parts = ax.violinplot(
                Y,
                [X],
                showmeans=False,
                showmedians=False,
                showextrema=False
            )
            for pc in parts['bodies']:
                pc.set_facecolor(C)
                pc.set_edgecolor('black')
                pc.set_alpha(1.0)

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        # ax.set_ylabel('density', fontsize=16)
        ax.set_ylabel(names[name]['xtitle'], fontsize=16)
        ax.set_ylim(names[name]['xlim'])
        ax.set_xlim(0, 5)
        ax.set_xticks([1, 2, 3, 4])
        ax.set_xticklabels(['a', 'b', 'c', 'd'])

        # ax.legend(fontsize=16)

        fig.tight_layout()
        fig.savefig(
            f'all_isomerdist_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


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
