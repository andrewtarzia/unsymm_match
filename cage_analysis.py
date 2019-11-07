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


def isomer_plot(dictionary, file_name, ytitle, ylim, horiz=None):
    """
    Get xTB energies of all four isomers.

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
            s=100
        )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xticks([X_positions[i] for i in X_positions])
    ax.set_xticklabels(list(X_positions.keys()))
    ax.set_xlim(0, 10)
    ax.set_ylim(ylim)
    if horiz is not None:
        for i in horiz:
            ax.axhline(y=i[0], c=i[1][0], lw=2)
    fig.tight_layout()
    fig.savefig(
        file_name,
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def get_energy(name, cage):
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
        # Extract energy.
        xtb_energy = stk.XTBEnergy(
            xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
            output_dir=f'cage_ey_{name}',
            num_cores=6,
            charge=4,
            num_unpaired_electrons=0,
            electronic_temperature=300,
            unlimited_memory=True
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
        expts = {
            'li1_lk2_li5': (10, 'k'),
            'li2_lk2_li6': (12, 'r'),
            'li1_lk2_li4': (14, 'k'),
            'li4_lk2_li5': (16, 'r')
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
        print(name_)
        energies[iso] = get_energy(name=name_, cage=cages[iso])

    min_energy = min(energies.values())
    energies = {
        i: energies[i]-min_energy for i in energies
    }

    isomer_plot(
        dictionary=energies,
        file_name=f'{name}_energies_plot.pdf',
        ytitle=r'relative energies [kJ/mol]',
        ylim=(-5, 300),
        horiz=experimental_lines()
    )
    return energies
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
            r'sum of plane deviation [$\mathrm{\AA}$]',
            (0, 2)
        )
    }

    for iso in cages:
        name_ = f'{name}_{iso}'
        print(name_)
        cage = cages[iso]
        results = atools.get_square_planar_distortion(
            mol=cage,
            metal=46,
            bonder=7
        )
        print(results)
        for measure in results:
            print(measure)
            print(results[measure])
            m_distortions[measure][0][iso] = np.mean(results[measure])
            print(m_distortions[measure][0][iso])

    for i in m_distortions:
        isomer_plot(
            dictionary=m_distortions[i][0],
            file_name=f'{name}_{m_distortions[i][1]}_plot.pdf',
            ytitle=m_distortions[i][2],
            ylim=m_distortions[i][3]
        )
    return m_distortions
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
        print(energies, energy_sep, energy_cutoff)
        if energy_sep < energy_cutoff:
            return False
    else:
        return False

    return True
