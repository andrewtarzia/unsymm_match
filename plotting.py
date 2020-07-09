#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting cage properties.

Author: Andrew Tarzia

Date Created: 29 May 2020

"""

import matplotlib.pyplot as plt
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
            'xtitle': r'min. $q_{\mathrm{sp}}$',
        }
    }

    for name, data in zip(names, [plane_devs, sqpl_ops]):
        fig, ax = plt.subplots(figsize=(8, 5))

        c_passed = atools.colors_i_like()[4]
        c_failed = atools.colors_i_like()[3]
        c_negative = atools.colors_i_like()[2]

        for i, lig in enumerate(ligands):
            if lig in cages_cis_wins:
                c = c_passed
            elif lig in cages_not_wins:
                c = c_failed

            if lig in experiments:
                m = 'X'
            else:
                m = 'o'

            if energy_preferences[i] < 0:
                c = c_negative

            ax.scatter(
                data[i],
                energy_preferences[i],
                c=c,
                edgecolors='k',
                marker=m,
                alpha=1,
                s=70
            )
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(names[name]['xtitle'], fontsize=16)
        ax.set_ylabel('energy preference [kJ/mol]', fontsize=16)
        ax.set_xlim(names[name]['xlim'])
        ax.set_ylim(-80, 80)

        ax.axhline(y=6.0, c='k', alpha=0.2)

        ax.scatter(
            -100,
            0,
            c=c_passed,
            edgecolors='k',
            marker='o',
            alpha=1,
            s=70,
            label='passed'
        )
        ax.scatter(
            -100,
            0,
            c=c_failed,
            edgecolors='k',
            marker='o',
            alpha=1,
            s=70,
            label='failed'
        )
        ax.scatter(
            -100,
            0,
            c=c_negative,
            edgecolors='k',
            marker='o',
            alpha=1,
            s=70,
            label='cis not preferred'
        )
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
    ax.set_xlim(0, 85)

    ax.bar(
        x_passed,
        y_passed,
        color=c_passed,
        width=1,
        edgecolor='none',
        alpha=1,
        label='passed'
    )
    ax.bar(
        x_failed,
        y_failed,
        color=c_failed,
        width=1,
        edgecolor='none',
        alpha=1,
        label='failed'
    )
    ax.bar(
        x_experiments,
        y_experiments,
        color=c_experiments,
        width=1,
        edgecolor='none',
        alpha=1,
        label='experimental cases'
    )
    ax.legend(fontsize=16)
    ax.axhline(y=y_bar, c='k', alpha=0.2)

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
