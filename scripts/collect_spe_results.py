#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to collect output from Gaussian 16 DFT SPE calculations.

Author: Andrew Tarzia

Date Created: 1 Mar 2021

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
        mol_name = splits[1].lower()+'_'+splits[2][0].upper()
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
        '3D1': {}, '4D2': {}, '5B4': {}, '5D1': {}, '5D3': {},
        '4B1': {}, '4B3': {}, '4C1': {}, '4C3': {}, '5A1': {},
        '5A3': {},
    }

    for rt in runtypes:
        rtd = runtypes[rt]
        for ser in series_data:
            a_energy = rtd['energies'][f'{ser.lower()}_A']
            b_energy = rtd['energies'][f'{ser.lower()}_B']
            c_energy = rtd['energies'][f'{ser.lower()}_C']
            d_energy = rtd['energies'][f'{ser.lower()}_D']
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
        ax.set_ylim(None, None)
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig(
            f'spe_{ser}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def ey_table(runtypes, selected_ligands):

    for rt in runtypes:
        print('-----', rt)
        columns = {'isomer': ['a', 'b', 'c', 'd', 'delta']}
        for i in selected_ligands:
            if i in ['4C1', '4C3']:
                continue
            columns[i] = []
        rtd = runtypes[rt]
        for sl in selected_ligands:
            if sl in ['4C1', '4C3']:
                continue
            a_energy = rtd['energies'][f'{sl.lower()}_A']
            b_energy = rtd['energies'][f'{sl.lower()}_B']
            c_energy = rtd['energies'][f'{sl.lower()}_C']
            d_energy = rtd['energies'][f'{sl.lower()}_D']
            all_energies = [a_energy, b_energy, c_energy, d_energy]
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
        '4B1', '4B3', '5A1', '5A3', '5B4',
        '4C1', '4C3',
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

    print(runtypes)

    plot_all_comparisons(runtypes)
    print('-----')
    ey_table(
        runtypes=runtypes,
        selected_ligands=selected_ligands,
    )
    print('-----')


if __name__ == '__main__':
    main()
