#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write input for ORCA 4.2.1 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

import json
import glob


def get_spe_energy(dir, basename):

    out_file = f'{dir}/{basename}.out'

    with open(out_file, 'r') as f:
        lines = f.readlines()

    targ_str = 'FINAL SINGLE POINT ENERGY'
    for line in lines:
        if targ_str in line:
            energy = float(line.rstrip().split(' ')[-1])
            break

    return energy


def main():

    list_of_mols = sorted(glob.glob('*.mol'))
    grid_lists = [
        'Grid0', 'Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5', 'Grid6',
        'Grid7',
    ]

    energies = {}

    for moll in list_of_mols:
        prefix = moll.replace('.mol', '')
        energies[prefix] = {}
        for grid in grid_lists:
            calc_name = f'{prefix}_spe_{grid}'
            directory = f'dir_{prefix}'
            grid_int = int(grid.replace('Grid', ''))
            energies[prefix][grid_int] = {}
            energies[prefix][grid_int]['dir'] = directory
            energies[prefix][grid_int]['cname'] = calc_name
            energies[prefix][grid_int]['energy'] = get_spe_energy(
                dir=directory,
                basename=calc_name
            )

    with open('collated_energies.json', 'w') as f:
        json.dump(energies, f, indent=4)


if __name__ == '__main__':
    main()
