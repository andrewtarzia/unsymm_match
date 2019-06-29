#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot energy of isomers of a MOC with nonsymmetric linkers.

Author: Andrew Tarzia

Date Created: 22 Jun 2019

"""

import sys
import logging
import matplotlib.pyplot as plt


def energy_data():
    """
    Returns list of edge directions for all isomers.

    Based on stk.M2L4_Lantern() with two building blocks and includes
    the heteroleptic cages.

    Returns
    -------
    :class:`list`
        List of building block placements to define cage.

    """
    return {
        '1111': {
            'gas_xtb': -238.866118441457,
            'DMSO_xtb': -239.479582872951,
            'MeCN_xtb': -239.467189113410,
            'gas_pm6': 0,
            'DMSO_pm6': 0,
            'MeCN_pm6': 0,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
        },
        '1-111': {
            'gas_xtb': -238.875561746491,
            'DMSO_xtb': -239.488637300703,
            'MeCN_xtb': -239.476575081578,
            'gas_pm6': 0,
            'DMSO_pm6': 0,
            'MeCN_pm6': 0,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
        },
        '11-1-1': {
            'gas_xtb': -238.876692139787,
            'DMSO_xtb': -239.488138205011,
            'MeCN_xtb': -239.476614755594,
            'gas_pm6': 0,
            'DMSO_pm6': 0,
            'MeCN_pm6': 0,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
        },
        '1-11-1': {
            'gas_xtb': -238.882446230086,
            'DMSO_xtb': -239.495093054748,
            'MeCN_xtb': -239.483584358587,
            'gas_pm6': 0,
            'DMSO_pm6': 0,
            'MeCN_pm6': 0,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
        }
    }


def three_plots(solv, filename, energies):
    fig, ax = plt.subplots()
    X_positions = {
        'all': 2, '1-down': 4, 'cis': 6, 'trans': 8
    }
    cs = ['k', 'r', 'b']
    ms = ['o', 'X', 'P']
    for i, c, m in zip(solv, cs, ms):
        X_values = []
        Y_values = []
        for j, iso in enumerate(energies):
            dc = energies[iso]
            X_values.append(list(X_positions.values())[j])
            Y_values.append(dc[i] * 2625.50)
        new_Y_values = [i-min(Y_values) for i in Y_values]
        print(i, new_Y_values)
        ax.plot(
            X_values, new_Y_values, c=c, label=i, marker=m
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('isomer label', fontsize=16)
    ax.set_ylabel('rel. total energy [kJ/mol]', fontsize=16)
    ax.set_xlim(0, 10)
    ax.set_ylim(-10, 50)
    ax.set_xticks(list(X_positions.values()))
    ax.set_xticklabels(list(X_positions.keys()))
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: plot_direction_isomers.py
        """)
        sys.exit()

    energies = energy_data()

    solv = ['gas_xtb', 'DMSO_xtb', 'MeCN_xtb']
    filename = 'slide_3_isomer_energies_xtb.pdf'
    three_plots(solv, filename, energies)

    solv = ['gas_pm6', 'DMSO_pm6', 'MeCN_pm6']
    filename = 'slide_3_isomer_energies_pm6.pdf'
    three_plots(solv, filename, energies)

    solv = ['gas_dft', 'DMSO_dft', 'MeCN_dft']
    filename = 'slide_3_isomer_energies_dft.pdf'
    three_plots(solv, filename, energies)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
