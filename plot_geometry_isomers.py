#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot energy of isomers of a MOC with nonsymmetric linkers
based on geometry effects.

Author: Andrew Tarzia

Date Created: 16 Jul 2019

"""

import sys
import logging
import matplotlib.pyplot as plt


def energy_data():
    """
    Returns total free energy of all isomers from all methods.

    Returns
    -------
    :class:`list`
        List of building block placements to define cage.

    """
    return {
        'cage1': {
            '1111': {
                'gas - xtb': -263.490590460015,
                'DMSO - xtb': -264.105987052207,
            },
            '1112': {
                'gas - xtb': -263.495225539274,
                'DMSO - xtb': -264.111998704079,
            },
            '1122': {
                'gas - xtb': -263.497111563998,
                'DMSO - xtb': -264.113844235440,
            },
            '1212': {
                'gas - xtb': -263.496107977462,
                'DMSO - xtb': -264.112553678762,
            },
        },
        'cage2': {
            '1111': {
                'gas - xtb': -208.344281547373,
                'DMSO - xtb': -208.979647567300,
            },
            '1112': {
                'gas - xtb': -208.348232129357,
                'DMSO - xtb': -208.983862936801,
            },
            '1122': {
                'gas - xtb': -208.348822919629,
                'DMSO - xtb': -208.983973230052,
            },
            '1212': {
                'gas - xtb': -208.347359487125,
                'DMSO - xtb': -208.982314660766,
            },
        }
    }


def make_plot(filename, energies):
    fig, ax = plt.subplots()
    X_positions = {
        'all': (2, '1111'), '1-down': (4, '1112'),
        'cis': (6, '1122'), 'trans': (8, '1212')
    }
    cs = ['r', 'b']
    ms = ['o', 'X']
    solv = ['gas - xtb', 'DMSO - xtb']
    cages = ['cage1', 'cage2']
    for s, c in zip(solv, cs):
        for cage, m in zip(cages, ms):
            X_values = []
            Y_values = []
            for X in X_positions:
                X_values.append(X_positions[X][0])
                Y_values.append(
                    energies[cage][X_positions[X][1]][s] * 2625.50
                )
            new_Y_values = [i-min(Y_values) for i in Y_values]
            print(s, cage, new_Y_values)
            label = f'{cage} - {s}'
            ax.plot(
                X_values, new_Y_values, c=c, label=label, marker=m
            )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('isomer label', fontsize=16)
    ax.set_ylabel('rel. total free energy [kJ/mol]', fontsize=16)
    ax.set_xlim(0, 10)
    ax.set_ylim(-5, 30)
    ax.set_xticks([X_positions[i][0] for i in X_positions])
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

    filename = 'geometry_isomer_energies_xtb.pdf'
    make_plot(filename, energies)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
