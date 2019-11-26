#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for analysing cages.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

import matplotlib.pyplot as plt
from os.path import exists
import numpy as np

import stk
import atools


def plot_energetics_and_geom(
    ligands,
    experiments,
    cages_cis_wins,
    cages_not_wins,
    energy_preferences,
    plane_devs
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
            c = c_negative

        ax.scatter(
            plane_devs[i],
            energy_preferences[i],
            c=c,
            edgecolors='k',
            marker=m,
            alpha=1,
            s=70
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(
        r'avg. plane deviation [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_ylabel('energy preference [kJ/mol]', fontsize=16)
    ax.set_xlim(0, 1)
    ax.set_ylim(-100, 100)

    ax.axhline(y=7.5, c='k', alpha=0.2)

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
        'all_cages_pref_and_stable.pdf',
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
            print(lig, energy_preferences[i])
            c = c_negative

        if energy_preferences[i] > 50:
            print(lig, energy_preferences[i])

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

    ax.axhline(y=7.5, c='k', alpha=0.2)

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


def isomer_plot(dictionary, file_name, ytitle, ylim, horiz=None):
    """
    Generic plot of isomer properties.

    Parameters
    ----------

    Returns
    -------

    """

    X_positions = {'A': 2, 'B': 4, 'C': 6, 'D': 8}

    fig, ax = plt.subplots(figsize=(8, 5))
    for isomer in dictionary:
        Y = dictionary[isomer]
        ax.scatter(
            X_positions[isomer],
            Y,
            c=atools.colors_i_like()[6],
            edgecolors='k',
            marker='o',
            alpha=1,
            s=120
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
            ax.axhline(y=i, c=j, lw=2, alpha=0.2)
    fig.tight_layout()
    fig.savefig(
        file_name,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def get_energy(name, cage, solvent=None):
    """
    Get xTB energy of a cage.

    Parameters
    ----------

    Returns
    -------

    """

    energy_file = f'{name}_optc.ey'
    if exists(energy_file):
        # Read .ey file.
        with open(energy_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            energy = float(line.rstrip())
            break
    else:
        print(f'calculating energy of {name}')
        # Extract energy.
        if solvent is None:
            solvent_str = None
            solvent_grid = 'normal'
        else:
            solvent_str, solvent_grid = solvent
        xtb_energy = stk.XTBEnergy(
            xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
            output_dir=f'cage_ey_{name}',
            num_cores=6,
            charge=4,
            num_unpaired_electrons=0,
            electronic_temperature=300,
            unlimited_memory=True,
            solvent=solvent_str,
            solvent_grid=solvent_grid
        )
        energy = xtb_energy.get_energy(cage)
        # Save to .ey file.
        with open(energy_file, 'w') as f:
            f.write(f'{energy}\n')

    return energy*2625.5


def get_cage_energies(name, cages):
    """
    Get xTB energies of all four isomers.

    Parameters
    ----------

    Returns
    -------

    """

    def experimental_lines():
        # These lines are the energy difference between the cis and
        # next most stable isomer from xtb after RDKIT opt.
        expts = {
            'li1_lk2_li5': (17, 'k'),
            'li2_lk2_li6': (27, 'k'),
            'li1_lk2_li4': (8, 'k'),
            'li4_lk2_li5': (1.4, 'r')
        }
        lines = [
            expts['li1_lk2_li5'][0],
            expts['li2_lk2_li6'][0],
            expts['li1_lk2_li4'][0],
            expts['li4_lk2_li5'][0],
        ]
        styles = [
            expts['li1_lk2_li5'][1],
            expts['li2_lk2_li6'][1],
            expts['li1_lk2_li4'][1],
            expts['li4_lk2_li5'][1],
        ]
        return lines, styles

    energies = {'A': None, 'B': None, 'C': None, 'D': None}

    for iso in energies:
        name_ = f'{name}_{iso}'
        energies[iso] = get_energy(
            name=name_,
            cage=cages[iso],
            solvent=('DMSO', 'verytight')
        )

    min_energy = min(energies.values())
    energies = {
        i: energies[i]-min_energy for i in energies
    }

    isomer_plot(
        dictionary=energies,
        file_name=f'{name}_energies_plot.pdf',
        ytitle=r'relative energies [kJ/mol]',
        ylim=(-5, 50),
        horiz=experimental_lines()
    )
    return energies


def get_ligand_distortion(name, cages, NN_dists, bites_dist):
    """
    Get ligand distortions compared to free ligand.

    Parameters
    ----------

    Returns
    -------

    """

    l_distortions = {
        'bite_angle': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'deltabites',
            r'average($\Delta$bite angle) [degrees]',
            (0, 100)
        ),
        'NN_dist': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'deltaNNs',
            r'average($\Delta$NN distance) [$\mathrm{\AA}$]',
            (0, 10)
        ),
    }
    for iso in cages:
        name_ = f'{name}_{iso}'
        cage = cages[iso]
        results = atools.calculate_ligand_distortion(
            mol=cage,
            cage_name=name_,
            free_ligand_name=f'{name}_opt_chosen.mol',
            free_NN_dists=NN_dists,
            free_bite_dists=bites_dist,
        )
        NN_change, bite_change = results
        l_distortions['NN_dist'][0][iso] = NN_change
        l_distortions['bite_angle'][0][iso] = bite_change

    for i in l_distortions:
        isomer_plot(
            dictionary=l_distortions[i][0],
            file_name=f'{name}_{l_distortions[i][1]}_plot.pdf',
            ytitle=l_distortions[i][2],
            ylim=l_distortions[i][3]
        )
    return l_distortions


def get_metal_centre_distortion(name, cages):
    """
    Get metal centre distortions.

    Parameters
    ----------

    Returns
    -------

    """

    m_distortions = {
        'bond_lengths': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'bonds',
            r'mean bond length [$\mathrm{\AA}$]',
            (1, 3)
        ),
        'angles': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'angles',
            'mean N-Pd-N angle [degrees]',
            (0, 100)
        ),
        'torsions': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'torsions',
            'mean N-N-N-N torsion [degrees]',
            (0, 100)
        ),
        'plane_dev': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'planes',
            r'max. sum of plane deviation [$\mathrm{\AA}$]',
            (0, 2)
        )
    }

    for iso in cages:
        cage = cages[iso]
        results = atools.get_square_planar_distortion(
            mol=cage,
            metal=46,
            bonder=7
        )
        for measure in results:
            if measure == 'plane_dev':
                m_distortions[measure][0][iso] = max(results[measure])
            else:
                m_distortions[measure][0][iso] = np.mean(
                    results[measure]
                )

    for i in m_distortions:
        isomer_plot(
            dictionary=m_distortions[i][0],
            file_name=f'{name}_{m_distortions[i][1]}_plot.pdf',
            ytitle=m_distortions[i][2],
            ylim=m_distortions[i][3]
        )
    return m_distortions


def check_stability(l_distortions, m_distortions):
    """
    Check if cis isomer is stable based on distortions.

    Parameters
    ----------

    Returns
    -------

    """

    def fail_NN(value):
        return value > 1

    def fail_plane_dev(value):
        return value > 0.3

    checks = {
        'bond_lengths': None,
        'angles': None,
        'torsions': None,
        'plane_dev': fail_plane_dev,
        'bite_angle': None,
        'NN_dist': None
    }

    for i in l_distortions:
        # No test for this measure.
        if checks[i] is None:
            continue
        check = checks[i](l_distortions[i][0]['C'])
        print(
            f'> {i}:',
            round(l_distortions[i][0]['C'], 4),
            check
        )
        if checks[i](l_distortions[i][0]['C']):
            return False

    for i in m_distortions:
        # No test for this measure.
        if checks[i] is None:
            continue
        check = checks[i](m_distortions[i][0]['C'])
        print(
            f'> {i}:',
            round(m_distortions[i][0]['C'], 4),
            check
        )
        if check:
            return False

    return True


def check_preference(energies, energy_cutoff):
    """
    Check if cis isomer is preferred based on relative energetics.

    Parameters
    ----------

    Returns
    -------

    """

    if energies['C'] == 0:
        energy_sep = min([
            energies[i] for i in energies if energies[i] != 0
        ])
        print(
            '> energy test:',
            round(energies['C'], 4),
            round(energy_sep, 4),
            energy_sep < energy_cutoff
        )
        if energy_sep < energy_cutoff:
            return False, energy_sep
    else:
        print(
            '> energy test (failed):',
            round(energies['C'], 4),
        )
        return False, -energies['C']

    return True, energy_sep
