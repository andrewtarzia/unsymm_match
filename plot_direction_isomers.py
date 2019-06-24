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
            'gas': -238.866118441457,
            'DMSO': -239.479582872951,
            'MeCN': -239.467189113410
        },
        '1-111': {
            'gas': -238.875561746491,
            'DMSO': -239.488637300703,
            'MeCN': -239.476575081578
        },
        '11-1-1': {
            'gas': -238.876692139787,
            'DMSO': -239.488138205011,
            'MeCN': -239.476614755594
        },
        '1-11-1': {
            'gas': -238.882446230086,
            'DMSO': -239.495093054748,
            'MeCN': -239.483584358587
        }
    }


def main():
    if (not len(sys.argv) == 1):
        print("""
    Usage: plot_direction_isomers.py
        """)
        sys.exit()

    energies = energy_data()

    fig, ax = plt.subplots()

    X_positions = {
        'all': 2, '1-down': 4, 'cis': 6, 'trans': 8
    }

    solv = ['gas', 'DMSO', 'MeCN']
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
        'slide_3_isomer_energies.pdf',
        dpi=720,
        bbox_inches='tight'
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
