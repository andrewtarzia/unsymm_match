#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to collect output from Gaussian 16 DFT SPE calculations.

Author: Andrew Tarzia

Date Created: 25 Feb 2021

"""

import matplotlib.pyplot as plt
import glob
import pandas as pd


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
        if back_splitter2 is None:
            raise ValueError()
        energy = energy.split(back_splitter2)[0]
        return float(energy)  # a.u.


def collate_energies(file_list, functional):

    collated_energies = {}
    for file in file_list:
        splits = file.split('/')
        mol_name = splits[1]
        splits = mol_name.split('_')
        mol_name = splits[1]+'_'+splits[2][0]
        energy = get_g16_energy(
            filename=file,
            front_splitter=f'SCF Done:  E({functional}) =',
            back_splitter='A.U. after',
            back_splitter2=None,
        )
        collated_energies[mol_name] = energy

    return collated_energies


def plot_all_comparisons(runtypes):

    series_data = {
        '3D1': {}, '4D2': {}, '5B4': {}, '5D1': {},
    }

    for rt in runtypes:
        print(rt)
        rtd = runtypes[rt]
        print(rtd)
        for ser in series_data:
            a_energy = rtd['energies'][f'{ser}_A']
            b_energy = rtd['energies'][f'{ser}_B']
            c_energy = rtd['energies'][f'{ser}_C']
            d_energy = rtd['energies'][f'{ser}_D']
            series_data[ser][rt] = [
                a_energy, b_energy, c_energy, d_energy
            ]

    X_positions = {'A': 2, 'B': 4, 'C': 6, 'D': 8}
    for ser in series_data:
        sd = series_data[ser]
        fig, ax = plt.subplots(figsize=(8, 5))
        for rt in sd:
            x = [i for i in X_positions.values()]
            eys = sd[rt]
            y = [i-min(eys) for i in eys]
            y = [i*2625.5 for i in y]
            ax.plot(
                x,
                y,
                color=runtypes[rt]['c'],
                lw=3,
                linestyle='--',
                alpha=0.3
            )
            ax.scatter(
                x,
                y,
                c=runtypes[rt]['c'],
                edgecolors='none',
                marker='o',
                alpha=1,
                s=180,
                label=runtypes[rt]['name'],
            )

        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_ylabel('rel. energy [kJ/mol]', fontsize=16)
        ax.set_xticks([X_positions[i] for i in X_positions])
        ax.set_xticklabels(list(X_positions.keys()))
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 250)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig(
            f'spe_test_{ser}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def main():

    runtypes = {
        'xtb': {
            'name': 'GFN2-xTB',
            'functional': None,
            'c': 'k',
        },
        'pbesdd': {
            'name': 'PBE0/def2-SVP/SDD',
            'functional': 'RPBE1PBE',
            'c': 'gold',
        },
        'pbelan': {
            'name': 'PBE0/def2-SVP/LANL2DZ',
            'functional': 'RPBE1PBE',
            'c': 'orange',
        },
        'pbenoecp': {
            'name': 'PBE0/def2-SVP',
            'functional': 'RPBE1PBE',
            'c': 'skyblue',
        },
        'b3lyplan': {
            'name': 'B3LYP/def2-SVP/SDD',
            'functional': 'RB3LYP',
            'c': 'forestgreen',
        },
        'mo6': {
            'name': 'M062X/def2-TZVP',
            'functional': 'RM062X',
            'c': 'firebrick',
        },
    }

    print(runtypes)
    input('check runtypes')

    for dir in runtypes:
        if dir == 'xtb':
            full_data = pd.read_csv(
                '/data/atarzia/projects/unsymm/screening/'
                'production/all_cage_results.txt'
            )
            target_df = full_data[
                full_data['lig'].isin(['3D1', '4D2', '5B4', '5D1'])
            ]

            runtypes[dir]['energies'] = {}
            for idx, row in target_df.iterrows():
                a_energy = float(row['energy_A'])
                b_energy = float(row['energy_B'])
                c_energy = float(row['energy_C'])
                d_energy = float(row['energy_D'])
                a_name = row['lig'] + '_A'
                b_name = row['lig'] + '_B'
                c_name = row['lig'] + '_C'
                d_name = row['lig'] + '_D'
                runtypes[dir]['energies'][f'{a_name}'] = a_energy/2625.5
                runtypes[dir]['energies'][f'{b_name}'] = b_energy/2625.5
                runtypes[dir]['energies'][f'{c_name}'] = c_energy/2625.5
                runtypes[dir]['energies'][f'{d_name}'] = d_energy/2625.5

        else:
            files = sorted(glob.glob(f'{dir}/*log'))
            functional = runtypes[dir]['functional']
            runtypes[dir]['energies'] = collate_energies(
                files, functional
            )

    print(runtypes)

    plot_all_comparisons(runtypes)


if __name__ == '__main__':
    main()
