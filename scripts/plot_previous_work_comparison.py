#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write input for Gaussian 16 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

import pandas as pd
import glob
import matplotlib.pyplot as plt

from utilities import bar_chart, collate_energies, colors_i_like


def plot_bar_charts(runtypes, ligands):

    series_data = {i: {} for i in ligands}

    for rt in runtypes:
        rtd = runtypes[rt]
        for ser in series_data:
            a_energy = rtd['energies'][f'{ser.lower()}_A']
            b_energy = rtd['energies'][f'{ser.lower()}_B']
            c_energy = rtd['energies'][f'{ser.lower()}_C']
            d_energy = rtd['energies'][f'{ser.lower()}_D']
            energy_sep = c_energy - min([a_energy, b_energy, d_energy])
            energy_sep = -1 * energy_sep*2625.5
            series_data[ser][rt] = energy_sep

    bar_chart(
        energies=[series_data[s]['xtb'] for s in series_data],
        filename='previous_xtb.pdf',
        x_labels=[s for s in series_data],
        y_title=r'xTB $cis$ energy preference [kJmol$^{-1}$]',
        facecolor=colors_i_like()[4],
    )
    bar_chart(
        energies=[series_data[s]['pbenoecp'] for s in series_data],
        filename='previous_dft.pdf',
        x_labels=[s for s in series_data],
        y_title=r'DFT $cis$ energy preference [kJmol$^{-1}$]',
        facecolor=colors_i_like()[3],
    )


def plot_(runtypes, ligands):

    series_data = {i: {} for i in ligands}

    for rt in runtypes:
        rtd = runtypes[rt]
        for ser in series_data:
            a_energy = rtd['energies'][f'{ser.lower()}_A']
            b_energy = rtd['energies'][f'{ser.lower()}_B']
            c_energy = rtd['energies'][f'{ser.lower()}_C']
            d_energy = rtd['energies'][f'{ser.lower()}_D']
            energy_sep = c_energy - min([a_energy, b_energy, d_energy])
            energy_sep = -1 * energy_sep*2625.5
            series_data[ser][rt] = energy_sep

    x_data = [series_data[s]['pbenoecp'] for s in series_data]
    y_data = [series_data[s]['xtb'] for s in series_data]
    x_title = r'DFT $cis$ energy preference [kJmol$^{-1}$]'
    y_title = r'xTB $cis$ energy preference [kJmol$^{-1}$]'

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(
        x_data,
        y_data,
        color='gold',
        edgecolor='k',
        s=160,
        alpha=1.0,
    )
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(x_title, fontsize=16)
    ax.set_ylabel(y_title, fontsize=16)
    ax.set_xlim(-5, None)
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig('previous_scatter.pdf', dpi=720, bbox_inches='tight')
    plt.close()


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
        '3D1', '4D2', '5D1', '5D3',
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

    plot_bar_charts(runtypes, selected_ligands)
    plot_(runtypes, selected_ligands)


if __name__ == '__main__':
    main()
