#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to make DRAFT latex table of properties of the selected ligands.

Author: Andrew Tarzia

Date Created: 30 Jan 2021

"""

import pandas as pd


def latex_table_body(selected_ligands, data):
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
        columns['ligand'].append(row['lig'])
        columns['qscp_C'].append(round(float(row['sqpl_op_C']), 2))
        columns['pd_C'].append(round(float(row['plane_dev_C']), 2))
        energy_sep = min([
            float(row['energy_A']),
            float(row['energy_B']),
            float(row['energy_D']),
        ])
        columns['xtb-energy-sep'].append(round(energy_sep, 1))
        columns['dft-energy-sep'].append('-')
        columns['expt-outcome'].append('X')

    final_df = pd.DataFrame(columns)

    print(final_df.to_latex(index=False))


def main():

    selected_ligands = [
        '5D1', '4D2', '5D3', '3D1',
        '4B3', '4B1', '5B4',
        '5A3', '5A1',
        # '5C2', '4C2', '3C2',
    ]

    full_data = pd.read_csv('all_cage_results.txt')

    latex_table_body(
        selected_ligands=selected_ligands,
        data=full_data,
    )


if __name__ == "__main__":
    main()
