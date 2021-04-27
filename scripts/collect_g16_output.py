#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to collect output from Gaussian 16 DFT SPE calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

from os.path import join
import json


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
        energy = energy.split(back_splitter2)[0]
        return float(energy)  # a.u.


def collate_energies(file_list):

    collated_energies = {}
    for file in file_list:
        splits = file.split('_')
        dir = 'opt_'+splits[1]+'_'+splits[2][0]
        mol_name = splits[1]+'_'+splits[2][0]
        print(dir, mol_name)
        output_file = join(dir, file)
        print(output_file)
        energy = get_g16_energy(
            filename=output_file,
            front_splitter='SCF Done:  E(RPBE1PBE) =',
            back_splitter='A.U. after',
            back_splitter2='a.u. after',
        )
        print(energy)
        input()
        collated_energies[mol_name] = energy

    return collated_energies


def main():

    with_solvent_files = [
        'o3_3d1_b6.log', 'o3_3d1_c1.log',
        'o3_4b1_a7.log', 'o3_4b1_b7.log',
        'o3_4b1_c6.log', 'o3_4b1_d2.log', 'o3_4b3_a5.log',
        'o3_4b3_br4.log', 'o3_4b3_c5.log', 'o3_4b3_d2.log',
        'o3_4d2_a2.log', 'o3_4d2_b6.log', 'o3_4d2_c1.log',
        'o3_4d2_d1.log', 'o3_5a3_a4.log', 'o3_5a3_b2.log',
        'o3_5a3_c5.log', 'o3_5a3_d1.log',
    ]

    without_solvent_files = [
        'o2_3d1_a1.log', 'o2_3d1_b1.log', 'o2_3d1_c2.log',
        'o2_3d1_d1.log', 'o2_4b1_a1.log', 'o2_4b1_b1.log',
        'o2_4b1_c1.log', 'o2_4b1_d1.log', 'o2_4b3_ar1.log',
        'o2_4b3_b2.log', 'o2_4b3_c2.log', 'o2_4b3_d1.log',
        'o2_4c1_a2.log', 'o2_4c1_b2.log', 'o2_4c1_c1.log',
        'o2_4c1_d3.log', 'o2_4c3_b2.log', 'o2_4d2_a1.log',
        'o2_4d2_b1.log', 'o2_4d2_c2.log', 'o2_4d2_d1.log',
        'o2_5a1_a1.log', 'o2_5a1_b1.log', 'o2_5a1_c1.log',
        'o2_5a1_d1.log', 'o2_5a3_a1.log', 'o2_5a3_b1.log',
        'o2_5a3_c1.log', 'o2_5a3_d1.log', 'o2_5b4_a1.log',
        'o2_5b4_b1.log', 'o2_5b4_c1.log', 'o2_5b4_d1.log',
        'o2_5d1_a2.log', 'o2_5d1_b1.log', 'o2_5d1_c1.log',
        'o2_5d1_d2.log', 'o2_5d3_a1.log', 'o2_5d3_b1.log',
        'o2_5d3_c1.log', 'o2_5d3_dr1.log',
    ]

    with_solvent_energies = collate_energies(with_solvent_files)
    without_solvent_energies = collate_energies(without_solvent_files)

    print(with_solvent_energies, without_solvent_energies)

    with open('collated_energies.json', 'w') as f:
        json.dump(without_solvent_energies, f)
    with open('collated_solv_energies.json', 'w') as f:
        json.dump(with_solvent_energies, f)


if __name__ == '__main__':
    main()
