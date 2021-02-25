#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot bar chart of lewis2020 comparing xtb and dft.

Author: Andrew Tarzia

Date Created: 24 Feb 2021

"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import sys

from atools import colors_i_like


def plot_xtb_bar_chart(data, ligands):

    filename = 'previous_xtb_energies.pdf'
    y_title = r'xTB $cis$ energy preference [kJmol$^{-1}$]'

    energies = []
    for lig in ligands:
        lig_data = data[data['lig'] == lig]
        print(lig_data)
        # Can assume that C is lowest energy here.
        energies.append(min([
            float(lig_data['energy_A']),
            float(lig_data['energy_B']),
            float(lig_data['energy_D']),
        ]))

    bar_chart(
        energies=energies,
        filename=filename,
        x_labels=[i.lower() for i in ligands],
        y_title=y_title,
    )


def plot_dft_bar_chart(data, ligands):

    filename = 'previous_dft_energies.pdf'
    y_title = r'DFT $cis$ energy preference [kJmol$^{-1}$]'

    print(data)

    energies = []
    for lig in ligands:
        print(lig)
        li = lig.lower()
        # a_data = [data[i] for i in data if lig in i and '_A_' in i]
        # a_data = a_data if a_data != [] else [1E24]
        # b_data = [data[i] for i in data if lig in i and '_B_' in i]
        # b_data = b_data if b_data != [] else [1E24]
        # c_data = [data[i] for i in data if lig in i and '_C_' in i]
        # c_data = c_data if c_data != [] else [1E24]
        # d_data = [data[i] for i in data if lig in i and '_D_' in i]
        # d_data = d_data if d_data != [] else [1E24]
        a_data = data[f'{li}_a']
        b_data = data[f'{li}_b']
        c_data = data[f'{li}_c']
        d_data = data[f'{li}_d']
        if c_data == [1E24]:
            energies.append(0)
        else:
            energies.append(min([
                i - c_data
                for i in [a_data, b_data, d_data]
            ]))

    print([i*2625.5 for i in energies])

    bar_chart(
        energies=[i*2625.5 for i in energies],
        filename=filename,
        x_labels=[i.lower() for i in ligands],
        y_title=y_title,
    )


def plot_opposing_bar_chart(xtbdata, dftdata, ligands):

    filename = 'previous_oppos_energies.pdf'

    xtbenergies = []
    for lig in ligands:
        lig_data = xtbdata[xtbdata['lig'] == lig]
        print(lig_data)
        # Can assume that C is lowest energy here.
        xtbenergies.append(min([
            float(lig_data['energy_A']),
            float(lig_data['energy_B']),
            float(lig_data['energy_D']),
        ]))

    dftenergies = []
    for lig in ligands:
        li = lig.lower()
        a_data = dftdata[f'{li}_a']
        b_data = dftdata[f'{li}_b']
        c_data = dftdata[f'{li}_c']
        d_data = dftdata[f'{li}_d']
        print(lig, a_data, b_data, c_data, d_data)
        if c_data == 1E24:
            dftenergies.append(0)
        else:
            dftenergies.append(min([
                i - c_data
                for i in [a_data, b_data, d_data]
            ]))

    print([i*2625.5 for i in dftenergies])

    opposing_bar_chart(
        xtbenergies=xtbenergies,
        dftenergies=[i*2625.5 for i in dftenergies],
        filename=filename,
        x_labels=[i.lower() for i in ligands],
    )


def bar_chart(energies, filename, x_labels, y_title):

    xs = [1, 2, 3, 4]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(
        xs,
        energies,
        width=1,
        color='#900C3F',
        edgecolor='k',
        linewidth=2
    )
    ax.set_xticks(xs)
    ax.set_xticklabels(x_labels)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand', fontsize=16)
    ax.set_ylabel(y_title, fontsize=16)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
    plt.close()


def opposing_bar_chart(
    xtbenergies,
    dftenergies,
    filename,
    x_labels,
):

    c_expt = colors_i_like()[3]
    c_selec = colors_i_like()[4]

    xs = [1, 2, 3, 4]

    fig, ax = plt.subplots(nrows=2, sharex=True, figsize=(8, 10))
    ax[0].bar(
        xs,
        xtbenergies,
        width=1,
        color=c_selec,
        edgecolor='k',
        linewidth=1.5,
        label='GFN2-xTB'
    )
    ax[1].bar(
        xs,
        dftenergies,
        width=1,
        color=c_expt,
        edgecolor='k',
        linewidth=1.5,
        label='PBE0-D3/def2-SVP/SDD/PCM(DMSO) [CHKTHIS]'
    )

    ax[1].set_xticks(xs)
    ax[1].set_xticklabels(x_labels)
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[0].set_ylabel(
        r'xTB $cis$ energy preference [kJmol$^{-1}$]',
        fontsize=16
    )
    ax[0].set_xlim(0, 5)
    ax[0].set_ylim(None, None)
    # ax[0].legend(fontsize=16)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    ax[1].set_ylabel(
        r'DFT $cis$ energy preference [kJmol$^{-1}$]',
        fontsize=16
    )
    ax[1].set_xlim(0, 5)
    ax[1].set_ylim(None, None)
    # ax[1].legend(fontsize=16)

    ax[1].invert_yaxis()
    # ax[0].set(yticks=df_male_1['age'])
    # ax[0].xaxis.tick_right()

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
    plt.close()


def main():
    if (not len(sys.argv) == 3):
        print(
            """
            Usage:
                dft_file (str) - File with DFT energies.

                xtb_file (str) - File with xtb energies.

            """
        )
        sys.exit()
    else:
        dft_file = sys.argv[1]
        xtb_file = sys.argv[2]

    experimental_ligands = ['5D3', '5D1', '4D2', '3D1']

    screen_data = pd.read_csv(xtb_file)
    with open(dft_file, 'r') as f:
        dft_data = json.load(f)

    plot_xtb_bar_chart(screen_data, experimental_ligands)
    plot_dft_bar_chart(dft_data, experimental_ligands)
    plot_opposing_bar_chart(
        xtbdata=screen_data,
        dftdata=dft_data,
        ligands=experimental_ligands,
    )


if __name__ == '__main__':
    main()
