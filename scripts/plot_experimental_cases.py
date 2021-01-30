#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot properties of cages chosen for experimental testing.

Author: Andrew Tarzia

Date Created: 10 Aug 2020

"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from atools import colors_i_like


def get_c(b):
    if b:
        return colors_i_like()[4]
    else:
        return colors_i_like()[3]


def get_xtick(name):

    return name


def get_C_pd(row):

    return float(row['plane_dev_C'])


def get_C_op(row):

    return float(row['sqpl_op_C'])


def get_energy_sep(row):

    energies = [
        row['energy_A'],
        row['energy_B'],
        row['energy_C'],
        row['energy_D']
    ]

    assert min(energies) == 0
    second_min = min([i for i in energies if i != 0])
    return second_min


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: plot_experimental_cases.py
        """)
        sys.exit()

    raise Exception('This script is outdated.')

    experimental_cases = {
        '5D1': {'new': False, 'p': 1},
        '4D2': {'new': False, 'p': 2},
        '5D3': {'new': False, 'p': 3},
        '3D1': {'new': False, 'p': 4},
        '5A3': {'new': True, 'p': 5},
        '5A1': {'new': True, 'p': 6},
        '4B3': {'new': True, 'p': 7},
        '4B1': {'new': True, 'p': 8},
        '5C2': {'new': True, 'p': 9},
        '4C2': {'new': True, 'p': 10},
        '5D7': {'new': True, 'p': 11},
        '5B4': {'new': True, 'p': 12},
    }

    print(experimental_cases)

    data = pd.read_csv('all_cage_results.txt')
    data = data[data['lig'].isin(experimental_cases)]
    print(data)
    print(data.columns)

    plot_properties = {
        'energies': {
            'ylim': (0, 40),
            'ytitle': 'energy separation [kJmol$^{-1}$]',
            'data_f': get_energy_sep,
            'ycut': 6,
        },
        'pds': {
            'ylim': (0, 0.2),
            'ytitle': r'max. plane deviation [$\mathrm{\AA}$]',
            'data_f': get_C_pd,
            'ycut': 0.3,
        },
        'qsps': {
            'ylim': (0.9, 1),
            'ytitle': r'min. $q_{\mathrm{sp}}$',
            'data_f': get_C_op,
            'ycut': 0.95,
        }
    }

    fig, (ax0, ax1, ax2) = plt.subplots(
        3,
        1,
        figsize=(8, 10),
        sharex=True,
    )
    ax0_dict = plot_properties['energies']
    ax0_data = [
        (row['lig'], ax0_dict['data_f'](row))
        for i, row in data.iterrows()
    ]
    ax0.bar(
        x=[experimental_cases[i[0]]['p'] for i in ax0_data],
        height=[i[1] for i in ax0_data],
        width=0.9,
        color=[
            get_c(experimental_cases[i[0]]['new']) for i in ax0_data
        ],
        tick_label=[get_xtick(i[0]) for i in ax0_data],
    )
    ax0.axhline(y=ax0_dict['ycut'], c='k', lw=2)
    ax0.tick_params(axis='both', which='major', labelsize=16)
    ax0.set_ylabel(ax0_dict['ytitle'], fontsize=16)
    ax0.set_ylim(ax0_dict['ylim'])

    ax1_dict = plot_properties['pds']
    ax1_data = [
        (row['lig'], ax1_dict['data_f'](row))
        for i, row in data.iterrows()
    ]
    ax1.bar(
        x=[experimental_cases[i[0]]['p'] for i in ax1_data],
        height=[i[1] for i in ax1_data],
        width=0.9,
        color=[
            get_c(experimental_cases[i[0]]['new']) for i in ax1_data
        ],
        tick_label=[get_xtick(i[0]) for i in ax1_data],
    )
    ax1.axhline(y=ax1_dict['ycut'], c='k', lw=2)
    ax1.tick_params(axis='both', which='major', labelsize=16)
    ax1.set_ylabel(ax1_dict['ytitle'], fontsize=16)
    ax1.set_ylim(ax1_dict['ylim'])

    ax2_dict = plot_properties['qsps']
    ax2_data = [
        (row['lig'], ax2_dict['data_f'](row))
        for i, row in data.iterrows()
    ]
    ax2.bar(
        x=[experimental_cases[i[0]]['p'] for i in ax2_data],
        height=[i[1] for i in ax2_data],
        width=0.9,
        color=[
            get_c(experimental_cases[i[0]]['new']) for i in ax2_data
        ],
        tick_label=[get_xtick(i[0]) for i in ax2_data],
    )
    ax2.axhline(y=ax2_dict['ycut'], c='k', lw=2)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    ax2.set_xlabel('ligand', fontsize=16)
    ax2.set_ylabel(ax2_dict['ytitle'], fontsize=16)
    ax2.set_ylim(ax2_dict['ylim'])

    legend_elements = [
        Patch(
            facecolor=colors_i_like()[4],
            edgecolor='none',
            label='new'
        ),
        Patch(
            facecolor=colors_i_like()[3],
            edgecolor='none',
            label='published'
        ),
    ]
    plt.xticks(rotation=45)
    ax0.legend(handles=legend_elements, fontsize=16)
    fig.tight_layout()
    fig.savefig(
        f'experimental_tests.pdf',
        dpi=720,
        bbox_inches='tight'
    )


if __name__ == "__main__":
    main()
