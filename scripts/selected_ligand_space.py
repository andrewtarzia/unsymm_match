#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot properties of the selected ligands.

Author: Andrew Tarzia

Date Created: 10 Sep 2020

"""

import pandas as pd
import matplotlib.pyplot as plt
from atools import colors_i_like


def bar_figure(selected_ligands, experimental_ligands, data):

    c_expt = colors_i_like()[3]
    c_selec = colors_i_like()[4]

    x_expt = {}
    x_selec = {}

    for i, li in enumerate(selected_ligands):
        row = data[data['lig'] == li]
        print(row)
        print(row['energy_C'])
        energy_sep = min([
            float(row['energy_A']),
            float(row['energy_B']),
            float(row['energy_D']),
        ])
        C_qsp = float(row['sqpl_op_C'])
        short_name = li
        if li in experimental_ligands:
            x_expt[i] = (short_name, energy_sep, C_qsp)
        else:
            x_selec[i] = (short_name, energy_sep, C_qsp)

    x_ticks = [i for i in x_expt]+[i for i in x_selec]
    x_ticklabels = (
        [x_expt[i][0] for i in x_expt]+[x_selec[i][0] for i in x_selec]
    )
    print(x_ticks, x_ticklabels)

    width = 0.9
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8, 8))

    ax1.barh(
        y=[i for i in x_expt],
        width=[x_expt[i][1] for i in x_expt],
        height=width,
        facecolor=c_expt,
        edgecolor='none',
        alpha=1,
        label='published examples'
    )
    ax1.barh(
        y=[i for i in x_selec],
        width=[x_selec[i][1] for i in x_selec],
        height=width,
        facecolor=c_selec,
        edgecolor='none',
        alpha=1,
        label='selected cage ligands'
    )

    # Set number of ticks for x-axis
    ax1.tick_params(axis='both', which='major', labelsize=16)
    # ax.set_xlabel(names[name]['xtitle'], fontsize=16)
    # ax1.set_xticks(x_ticks)
    # ax1.set_xticklabels(x_ticklabels)
    ax1.set_xlabel('stability of C isomer [kJmol$^{-1}$]', fontsize=16)
    # ax.set_xlim(names[name]['xlim'])
    ax1.set_xlim(0, 30)

    ax1.axvline(x=6.0, c='k', alpha=0.6, lw=2)

    ax2.barh(
        y=[i for i in x_expt],
        height=width,
        width=[x_expt[i][2] for i in x_expt],
        facecolor=c_expt,
        edgecolor='none',
        alpha=1,
    )
    ax2.barh(
        y=[i for i in x_selec],
        height=width,
        width=[x_selec[i][2] for i in x_selec],
        facecolor=c_selec,
        edgecolor='none',
        alpha=1,
    )

    # Set number of ticks for x-axis
    ax2.tick_params(axis='both', which='major', labelsize=16)
    ax2.set_xticks([0.9, 1])
    ax2.set_xticklabels(['0.9', '1.0'])
    ax2.set_yticks(x_ticks)
    ax2.set_yticklabels(x_ticklabels, rotation=45)
    ax2.set_xlabel(r'$q_{\mathrm{sqp,min}}$', fontsize=16)
    ax2.set_xlim(0.9, 1)
    ax1.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        f'selected_ligands.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def single_bar_figure(selected_ligands, experimental_ligands, data):

    c_expt = colors_i_like()[3]
    c_selec = colors_i_like()[4]

    x_expt = {}
    x_selec = {}

    for i, li in enumerate(selected_ligands):
        row = data[data['lig'] == li]
        print(row, row['energy_C'])
        energy_sep = min([
            float(row['energy_A']),
            float(row['energy_B']),
            float(row['energy_D']),
        ])
        C_qsp = float(row['sqpl_op_C'])
        short_name = (
            li.replace(
                '_li', '-'
            ).replace('_lk', '-').replace('li', '')
        )
        if li in experimental_ligands:
            x_expt[i] = (short_name, energy_sep, C_qsp)
        else:
            x_selec[i] = (short_name, energy_sep, C_qsp)

    x_ticks = [i for i in x_expt]+[i for i in x_selec]
    x_ticklabels = (
        [x_expt[i][0] for i in x_expt]+[x_selec[i][0] for i in x_selec]
    )
    print(x_ticks, x_ticklabels)

    width = 0.9
    fig, ax = plt.subplots(figsize=(4, 8))

    ax.barh(
        y=[i for i in x_expt],
        width=[x_expt[i][1] for i in x_expt],
        height=width,
        facecolor=c_expt,
        edgecolor='k',
        alpha=1,
        label='published examples'
    )
    ax.barh(
        y=[i for i in x_selec],
        width=[x_selec[i][1] for i in x_selec],
        height=width,
        facecolor=c_selec,
        edgecolor='k',
        alpha=1,
        label='selected cage ligands'
    )

    ax.axvline(x=6.0, c='k', alpha=0.6, lw=2)

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.set_xlabel(names[name]['xtitle'], fontsize=16)
    # ax1.set_xticks(x_ticks)
    # ax1.set_xticklabels(x_ticklabels)
    ax.set_xlabel('stability of C isomer [kJmol$^{-1}$]', fontsize=16)
    # ax.set_xlim(names[name]['xlim'])
    ax.set_xlim(0, 30)
    # Set number of ticks for x-axis
    ax.set_yticks(x_ticks)
    ax.set_yticklabels(x_ticklabels)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        f'selected_ligands_energysep.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def main():

    experimental_ligands = ['5D1', '4D2', '5D3', '3D1']
    selected_ligands = [
        '5D1', '4D2', '5D3', '3D1',
        '4B3', '4B1', '5B4',
        '5A3', '5A1',
        # '5C2', '4C2', '3C2',
    ]

    full_data = pd.read_csv('all_cage_results.txt')
    print(full_data.columns)

    bar_figure(
        selected_ligands=selected_ligands,
        experimental_ligands=experimental_ligands,
        data=full_data,
    )
    single_bar_figure(
        selected_ligands=selected_ligands,
        experimental_ligands=experimental_ligands,
        data=full_data,
    )


if __name__ == "__main__":
    main()
