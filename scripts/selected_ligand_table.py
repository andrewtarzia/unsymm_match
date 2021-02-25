#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to make DRAFT latex table of properties of the selected ligands.

Author: Andrew Tarzia

Date Created: 30 Jan 2021

"""

import pandas as pd
import json


def latex_table_body(selected_ligands, data, dft_data):
    print(data.columns)
    target_df = data[data['lig'].isin(selected_ligands)]

    columns = {
        'ligand': [],
        'qscp_C': [],
        'pd_C': [],
        'xtb-energy-sep': [],
        'dft-energy-sep': [],
        'expt-outcome': [],
    }
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
        a_data = dft_data[f'{li}_a']
        b_data = dft_data[f'{li}_b']
        c_data = dft_data[f'{li}_c']
        d_data = dft_data[f'{li}_d']
        dft_energy_sep = min([
            i - c_data
            for i in [a_data, b_data, d_data]
        ])
        dft_energy_sep = dft_energy_sep*2625.5

        columns['dft-energy-sep'].append(round(dft_energy_sep, 1))
        columns['expt-outcome'].append('X')

    final_df = pd.DataFrame(columns)

    print(final_df.to_latex(index=False))


def dft_table(selected_ligands, data, dft_data):

    target_df = data[data['lig'].isin(selected_ligands)]

    columns = {'isomer': ['a', 'b', 'c', 'd', 'delta']}
    for i in selected_ligands:
        columns[i] = []

    for idx, row in target_df.iterrows():
        li = row['lig'].lower()
        a_data = float(dft_data[f'{li}_a'])
        b_data = float(dft_data[f'{li}_b'])
        c_data = float(dft_data[f'{li}_c'])
        d_data = float(dft_data[f'{li}_d'])

        all_energies = [a_data, b_data, c_data, d_data]
        a_energy = a_data - min(all_energies)
        b_energy = b_data - min(all_energies)
        c_energy = c_data - min(all_energies)
        d_energy = d_data - min(all_energies)

        energy_sep = c_energy - min([a_energy, b_energy, d_energy])
        energy_sep = energy_sep*2625.5

        a_energy = a_energy*2625.5
        b_energy = b_energy*2625.5
        c_energy = c_energy*2625.5
        d_energy = d_energy*2625.5
        print(a_energy, b_energy, c_energy, d_energy)
        print(energy_sep)
        input()
        columns[row['lig']] = [
            round(a_energy, 1),
            round(b_energy, 1),
            round(c_energy, 1),
            round(d_energy, 1),
            round(energy_sep, 1),
        ]

    final_df = pd.DataFrame(columns)

    print(final_df.to_latex(index=False))


def xtb_table(selected_ligands, data):
    print(data.columns)

    target_df = data[data['lig'].isin(selected_ligands)]

    columns = {'isomer': ['a', 'b', 'c', 'd', 'delta']}
    for i in selected_ligands:
        columns[i] = []

    for idx, row in target_df.iterrows():
        a_energy = float(row['energy_A'])
        b_energy = float(row['energy_B'])
        c_energy = float(row['energy_C'])
        d_energy = float(row['energy_D'])
        print(a_energy, b_energy, c_energy, d_energy)
        energy_sep = c_energy - min([a_energy, b_energy, d_energy])
        print(energy_sep)
        columns[row['lig']] = [
            round(a_energy, 1),
            round(b_energy, 1),
            round(c_energy, 1),
            round(d_energy, 1),
            round(energy_sep, 1),
        ]

    final_df = pd.DataFrame(columns)

    print(final_df.to_latex(index=False))


def main():

    selected_ligands = [
        '3D1', '4D2', '5D1', '5D3',
        '4B1', '4B3', '5A1', '5A3', '5B4',
        # '5C2', '4C2', '3C2',
    ]

    full_data = pd.read_csv('all_cage_results.txt')
    with open('dft_validation/collated_energies.json', 'r') as f:
        dft_data = json.load(f)

    latex_table_body(
        selected_ligands=selected_ligands,
        data=full_data,
        dft_data=dft_data,
    )
    print('-----')
    xtb_table(
        selected_ligands=selected_ligands,
        data=full_data,
    )
    print('-----')
    dft_table(
        selected_ligands=selected_ligands,
        data=full_data,
        dft_data=dft_data,
    )


if __name__ == "__main__":
    main()
