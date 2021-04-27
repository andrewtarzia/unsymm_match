#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to extract energies from ORCA 4.2.1 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

import json
from os.path import exists
import glob
import matplotlib.pyplot as plt

from atools import replace


def get_spe_energy(dir, outfile):

    out_file = f'{dir}/{outfile}'

    with open(out_file, 'r') as f:
        lines = f.readlines()

    targ_str = 'FINAL SINGLE POINT ENERGY'
    for line in lines:
        if targ_str in line:
            energy = float(line.rstrip().split(' ')[-1])
            break

    return energy


def get_isomer(isomers, string):

    for i in isomers:
        if i in string:
            return i


def plot_isomers(edict):

    isomers = ['A', 'B', 'C', 'D']
    seen_prefix = set()
    isomer_energies = {}
    for i in edict:
        pref = replace(
            string=i,
            substitutions={
                '_optc': '',
                '_A': '',
                '_B': '',
                '_C': '',
                '_D': '',
            }
        )
        isomer = get_isomer(isomers=isomers, string=i)
        if pref in seen_prefix:
            isomer_energies[pref][isomer] = (
                float(edict[i]['energy'])*2625.5
            )
        else:
            isomer_energies[pref] = {}
            isomer_energies[pref][isomer] = (
                float(edict[i]['energy'])*2625.5
            )
            seen_prefix.add(pref)

    print(isomer_energies)
    X_positions = {'A': 2, 'B': 4, 'C': 6, 'D': 8}
    for cage in isomer_energies:
        print(cage)
        cage_dict = isomer_energies[cage]
        print(cage_dict)
        min_energy = min([
            cage_dict[iso] for iso in cage_dict
        ])
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(
            [X_positions[iso] for iso in cage_dict],
            [cage_dict[iso]-min_energy for iso in cage_dict],
            c='k',
            marker='o',
        )

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('isomer', fontsize=16)
        ax.set_ylabel('rel. total energy [kJmol$^{-1}$]', fontsize=16)
        ax.set_xticks([X_positions[i] for i in X_positions])
        ax.set_xticklabels(list(X_positions.keys()))
        ax.set_ylim(-1, 500)
        fig.tight_layout()
        fig.savefig(
            f'{cage}_spe_energy.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def main():

    list_of_mols = sorted(glob.glob('*.mol'))

    energies = {}
    for moll in list_of_mols:
        prefix = moll.replace('.mol', '')

        # for grid in grid_lists:
        calc_name = f'{prefix}'
        spe_directory = f'spe_{prefix}'
        print(f'> extracting {calc_name}.....')
        infile = f's_{calc_name}.in'
        # Should check here for different naming based on optimisation
        # steps.
        outfile = infile.replace('.in', '.out')
        if not exists(f'{spe_directory}/{outfile}'):
            raise FileNotFoundError(
                f'{spe_directory}/{outfile} not found. '
                'Run calculation.'
            )

        energies[prefix] = {}
        energies[prefix]['dir'] = spe_directory
        energies[prefix]['cname'] = calc_name
        energies[prefix]['energy'] = get_spe_energy(
            dir=spe_directory,
            outfile=outfile
        )

    with open('collated_spe_energies.json', 'w') as f:
        json.dump(energies, f, indent=4)

    plot_isomers(edict=energies)


if __name__ == '__main__':
    main()
