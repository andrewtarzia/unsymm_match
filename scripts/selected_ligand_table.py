#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to make DRAFT latex table of properties of the selected ligands.

Author: Andrew Tarzia

Date Created: 30 Jan 2021

"""

import pandas as pd
import glob

from utilities import collate_energies


def latex_table_body(full_data, runtypes, selected_ligands):
    target_df = full_data[full_data['lig'].isin(selected_ligands)]

    columns = {
        'ligand': [],
        'qscp_C': [],
        'pd_C': [],
        'xtb-energy-sep': [],
        'dft-energy-sep': [],
        'expt-outcome': [],
    }

    dft_data = runtypes['pbenoecp']['energies']
    print(dft_data)

    for idx, row in target_df.iterrows():
        columns['ligand'].append(row['lig'].lower())
        columns['qscp_C'].append(round(float(row['sqpl_op_C']), 3))
        columns['pd_C'].append(round(float(row['plane_dev_C']), 3))
        energy_sep = min([
            float(row['energy_A']),
            float(row['energy_B']),
            float(row['energy_D']),
        ])
        columns['xtb-energy-sep'].append(round(energy_sep, 1))
        li = row['lig'].lower()
        a_data = dft_data[f'{li}_A']
        b_data = dft_data[f'{li}_B']
        c_data = dft_data[f'{li}_C']
        d_data = dft_data[f'{li}_D']
        dft_energy_sep = min([
            i - c_data
            for i in [a_data, b_data, d_data]
        ])
        dft_energy_sep = dft_energy_sep*2625.5

        columns['dft-energy-sep'].append(round(dft_energy_sep, 1))
        columns['expt-outcome'].append('X')

    final_df = pd.DataFrame(columns)

    print(final_df.to_latex(index=False))


def main():

    runtypes = {
        'xtb': {
            'name': 'GFN2-xTB',
            'functional': None,
            'c': 'k',
        },
        'pbenoecp': {
            'name': 'PBE0/def2-SVP',
            'functional': 'RPBE1PBE',
            'c': 'skyblue',
        },
    }

    selected_ligands = [
        # '3D1', '4D2', '5D1', '5D3',
        '4B1', '4B3', '5A1', '5A3', '5B4',
    ]

    for dir in runtypes:
        if dir == 'xtb':
            full_data = pd.read_csv(
                '/data/atarzia/projects/unsymm/screening/'
                'production/all_cage_results.txt'
            )
            target_df = full_data[
                full_data['lig'].isin(selected_ligands)
            ]

            runtypes[dir]['energies'] = {}
            for idx, row in target_df.iterrows():
                a_energy = float(row['energy_A'])
                b_energy = float(row['energy_B'])
                c_energy = float(row['energy_C'])
                d_energy = float(row['energy_D'])
                a_name = row['lig'].lower() + '_A'
                b_name = row['lig'].lower() + '_B'
                c_name = row['lig'].lower() + '_C'
                d_name = row['lig'].lower() + '_D'
                runtypes[dir]['energies'][f'{a_name}'] = a_energy/2625.5
                runtypes[dir]['energies'][f'{b_name}'] = b_energy/2625.5
                runtypes[dir]['energies'][f'{c_name}'] = c_energy/2625.5
                runtypes[dir]['energies'][f'{d_name}'] = d_energy/2625.5

        else:
            files = sorted(glob.glob('*/*log'))
            functional = runtypes[dir]['functional']
            runtypes[dir]['energies'] = collate_energies(
                files, functional
            )

    latex_table_body(
        full_data=full_data,
        runtypes=runtypes,
        selected_ligands=selected_ligands,
    )


if __name__ == "__main__":
    main()
