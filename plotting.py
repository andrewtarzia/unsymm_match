#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting cage properties.

Author: Andrew Tarzia

Date Created: 29 May 2020

"""

import matplotlib.pyplot as plt
from pandas import read_csv
import numpy as np
from os.path import exists
from os import mkdir
from rdkit.Chem import AllChem as rdkit

import atools


def plot_energetics_and_geom(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    energy_preferences,
    plane_devs,
    sqpl_ops
):
    """
    Plot energy preference of all cages.

    Parameters
    ----------

    Returns
    -------

    """

    names = {
        'plane_dev': {
            'xlim': (0, 1),
            'xtitle': r'max. plane deviation [$\mathrm{\AA}$]',
        },
        'sqpl': {
            'xlim': (0.6, 1),
            'xtitle': r'$q_{\mathrm{sqp,min}}$',
        }
    }

    for name, data in zip(names, [plane_devs, sqpl_ops]):
        fig, ax = plt.subplots(figsize=(8, 5))

        c_passed = atools.colors_i_like()[4]
        # c_failed = atools.colors_i_like()[3]
        c_negative = atools.colors_i_like()[3]
        c_experiments = atools.colors_i_like()[0]
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
        ax.set_ylabel('stability of C isomer [kJ/mol]', fontsize=16)
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

    Parameters
    ----------

    Returns
    -------

    """

    fig, ax = plt.subplots(figsize=(8, 5))

    c_passed = atools.colors_i_like()[4]
    c_failed = atools.colors_i_like()[3]
    c_experiments = atools.colors_i_like()[0]

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


def draw_molecules(
    ligands,
    experiments,
    energy_preferences,
    plane_devs,
    sqpl_ops,
):
    """
    Draw molecules as grids with scores.

    Parameters
    ----------

    Returns
    -------

    """

    mol_list = []
    name_list = []
    count = 0
    for i, lig in enumerate(ligands):
        MOL = rdkit.MolFromMolFile(
            f'{lig}_opt.mol'
        )
        MOL.RemoveAllConformers()
        mol_list.append(MOL)
        name_list.append(
            f'{lig}: OP={round(sqpl_ops[i], 2)}'
            f',PD={round(plane_devs[i], 2)}\n'
            f'E={round(energy_preferences[i], 2)}'
        )
        count += 1

    # Sort by energy preferences.
    mol_list = [
        x for _, x in sorted(zip(energy_preferences, mol_list))
    ]
    name_list = [
        x for _, x in sorted(zip(energy_preferences, name_list))
    ]

    if not exists('molecules_scoring'):
        mkdir('molecules_scoring')

    # Save figure of desired molecules.
    atools.mol_list2grid(
        molecules=mol_list,
        names=name_list,
        filename='molecules_scoring/molecules_scoring',
        mol_per_row=4,
        maxrows=3,
        subImgSize=(300, 300)
    )


def isomer_plot(dictionary, file_name, ytitle, ylim, horiz=None):
    """
    Generic plot of isomer properties.

    Parameters
    ----------

    Returns
    -------

    """

    col = 3

    X_positions = {'A': 2, 'B': 4, 'C': 6, 'D': 8}

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        list(X_positions.values()),
        list(dictionary.values()),
        color=atools.colors_i_like()[col],
        lw=3,
    )
    for isomer in dictionary:
        Y = dictionary[isomer]
        ax.scatter(
            X_positions[isomer],
            Y,
            c=atools.colors_i_like()[col],
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

    Parameters
    ----------

    Returns
    -------

    """

    data = read_csv('all_cage_results.txt')

    names = {
        'plane_dev': {
            'xlim': (0, 1),
            'xtitle': r'max. plane deviation [$\mathrm{\AA}$]',
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
            'xtitle': r'relative isomer energy [kJ/mol]',
            'atitle': 'energy_A',
            'btitle': 'energy_B',
            'ctitle': 'energy_C',
            'dtitle': 'energy_D',
            'width': 5,
        },
    }

    for name in names:
        fig, ax = plt.subplots(figsize=(5, 5))
        X_bins = np.arange(
            names[name]['xlim'][0],
            names[name]['xlim'][1],
            names[name]['width'],
        )

        c_a = atools.colors_i_like()[4]
        c_b = atools.colors_i_like()[3]
        c_c = atools.colors_i_like()[0]
        c_d = atools.colors_i_like()[2]
        x_a = list(data[names[name]['atitle']])
        x_b = list(data[names[name]['btitle']])
        x_c = list(data[names[name]['ctitle']])
        x_d = list(data[names[name]['dtitle']])
        labels = ['A', 'B', 'C', 'D']
        XS = [1, 2, 3, 4]
        YS = [x_a, x_b, x_c, x_d]
        CS = [c_a, c_b, c_c, c_d]

        for X, label, Y, C in zip(XS, labels, YS, CS):
            hist, bin_edges = np.histogram(
                a=Y,
                bins=X_bins,
                density=True
            )
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
        ax.set_xticklabels(['A', 'B', 'C', 'D'])

        # ax.legend(fontsize=16)

        fig.tight_layout()
        fig.savefig(
            f'all_isomerdist_{name}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()
