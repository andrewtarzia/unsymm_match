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


def plot_energetics(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    energy_preferences
):
    """
    Plot energy preference of all cages.

    Parameters
    ----------

    Returns
    -------

    """

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
            print('negative:', lig, energy_preferences[i])
            c = c_negative

        if energy_preferences[i] > 50:
            print('really large:', lig, energy_preferences[i])

        ax.scatter(
            i+1,
            energy_preferences[i],
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1,
            s=70
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('cage', fontsize=16)
    ax.set_ylabel('energy preference [kJ/mol]', fontsize=16)
    ax.set_xlim(0, 62)

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
        'all_cages_energy_preferences.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_lses(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    lses
):
    """
    Plot sum of ligand strain energy of all cages.

    Parameters
    ----------

    Returns
    -------

    """

    c_passed = atools.colors_i_like()[4]
    c_failed = atools.colors_i_like()[3]

    fig, ax = plt.subplots(figsize=(8, 5))

    for i, lig in enumerate(ligands):
        if lig in cages_cis_wins:
            c = c_passed
        elif lig in cages_not_wins:
            c = c_failed

        if lig in experiments:
            m = 'X'
        else:
            m = 'o'

        ax.scatter(
            i+1,
            lses[i],
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1,
            s=70
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('cage', fontsize=16)
    ax.set_ylabel(
        r'sum(ligand strain energy) [kJmol$^{-1}$]',
        fontsize=16
    )
    ax.set_xlim(0, 62)

    print('fill in!')
    ax.axhline(y=0.3, c='k', alpha=0.2)

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
    ax.legend(fontsize=16)
    # ax.set_ylim(-1000, 1000)
    # if horiz is not None:
    #     for i, j in zip(*horiz):
    #         ax.axhline(y=i, c=j, lw=2, alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        'all_cages_LSES.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_plane_devs(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    plane_devs
):
    """
    Plot plane deviation of all cages.

    Parameters
    ----------

    Returns
    -------

    """

    c_passed = atools.colors_i_like()[4]
    c_failed = atools.colors_i_like()[3]

    fig, ax = plt.subplots(figsize=(8, 5))

    for i, lig in enumerate(ligands):
        if lig in cages_cis_wins:
            c = c_passed
        elif lig in cages_not_wins:
            c = c_failed

        if lig in experiments:
            m = 'X'
        else:
            m = 'o'

        ax.scatter(
            i+1,
            plane_devs[i],
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1,
            s=70
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('cage', fontsize=16)
    ax.set_ylabel(
        r'max. plane deviation [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_xlim(0, 62)

    ax.axhline(y=0.3, c='k', alpha=0.2)

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
    ax.legend(fontsize=16)
    # ax.set_ylim(-1000, 1000)
    # if horiz is not None:
    #     for i, j in zip(*horiz):
    #         ax.axhline(y=i, c=j, lw=2, alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        'all_cages_planedevs.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_sqpl_ops(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    sqpl_ops
):
    """
    Plot plane deviation of all cages.

    Parameters
    ----------

    Returns
    -------

    """

    c_passed = atools.colors_i_like()[4]
    c_failed = atools.colors_i_like()[3]

    fig, ax = plt.subplots(figsize=(8, 5))

    for i, lig in enumerate(ligands):
        if lig in cages_cis_wins:
            c = c_passed
        elif lig in cages_not_wins:
            c = c_failed

        if lig in experiments:
            m = 'X'
        else:
            m = 'o'

        ax.scatter(
            i+1,
            sqpl_ops[i],
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1,
            s=70
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel('cage', fontsize=16)
    ax.set_ylabel(
        r'min. $q_{\mathrm{sp}}$',
        fontsize=16
    )
    ax.set_xlim(0, 62)
    ax.set_ylim(0, 1)

    print('fill in!')
    ax.axhline(y=0.3, c='k', alpha=0.2)

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
    ax.legend(fontsize=16)
    # ax.set_ylim(-1000, 1000)
    # if horiz is not None:
    #     for i, j in zip(*horiz):
    #         ax.axhline(y=i, c=j, lw=2, alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        'all_cages_sqplop.pdf',
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
