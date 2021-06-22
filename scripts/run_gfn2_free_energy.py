#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run GFN2-xTB free energy calculations.

Author: Andrew Tarzia

Date Created: 21 Jun 2021

"""

import pandas as pd
import json
import os
import stk
import stko


def ey_table(all_results, selected_ligands):

    print(all_results)
    columns = {'isomer': ['a', 'b', 'c', 'd', 'delta']}
    for sl in selected_ligands:
        if sl in ['4C1', '4C3']:
            continue
        a_energy = all_results[f'{sl.lower()}_a']['fe']
        b_energy = all_results[f'{sl.lower()}_b']['fe']
        c_energy = all_results[f'{sl.lower()}_c']['fe']
        d_energy = all_results[f'{sl.lower()}_d']['fe']
        all_energies = [a_energy, b_energy, c_energy, d_energy]
        print(sl, all_energies)
        energy_sep = c_energy - min([a_energy, b_energy, d_energy])
        energy_sep = energy_sep*2625.5
        a_energy = a_energy - min(all_energies)
        a_energy = a_energy * 2625.5
        b_energy = b_energy - min(all_energies)
        b_energy = b_energy * 2625.5
        c_energy = c_energy - min(all_energies)
        c_energy = c_energy * 2625.5
        d_energy = d_energy - min(all_energies)
        d_energy = d_energy * 2625.5
        columns[sl] = [
            round(a_energy, 1),
            round(b_energy, 1),
            round(c_energy, 1),
            round(d_energy, 1),
            round(energy_sep, 1),
        ]
    final_df = pd.DataFrame(columns)

    print(final_df.to_latex(index=False))
    print('-----')


def main():
    num_proc = 6
    optc_extension = '_optc.mol'
    isomers = ['A', 'B', 'C', 'D']
    selected_ligands = [
        '5D1', '4D2', '5D3', '3D1',
        '4B3', '4B1', '5B4',
        '5A3', '5A1',
        '4C1', '4C3',
        # '5C2', '4C2', '3C2',
    ]

    # full_data = pd.read_csv('../all_cage_results.txt')
    # df_to_run = full_data[full_data['lig'].isin(selected_ligands)]

    all_results = {}
    for lig in selected_ligands:
        for isomer in isomers:
            struct_file = f'../{lig}_{isomer}{optc_extension}'
            prefix = f'{lig.lower()}_{isomer.lower()}'
            struct = stk.BuildingBlock.init_from_file(struct_file)
            opt_directory = f'fe_{prefix}'
            energy_file = f'{prefix}.fey'
            if os.path.exists(energy_file):
                # Read .ey file.
                with open(energy_file, 'r') as f:
                    res_dict = json.load(f)
            else:

                xtb = stko.XTBEnergy(
                    xtb_path=(
                        '/home/atarzia/software/xtb-190806/bin/xtb'
                    ),
                    output_dir=opt_directory,
                    unlimited_memory=True,
                    num_cores=num_proc,
                    calculate_free_energy=True,
                    charge=4,
                    solvent='DMSO',
                    solvent_grid='verytight',
                )

                xtb_results = xtb.get_results(struct)
                total_free_energy = xtb_results.get_total_free_energy()
                res_dict = {
                    'fe': total_free_energy[0],
                    'unit': total_free_energy[1],
                }
                # Save to .ey file.
                with open(energy_file, 'w') as f:
                    json.dump(res_dict, f)

            all_results[prefix] = res_dict

    ey_table(all_results, selected_ligands)


if __name__ == '__main__':
    main()
