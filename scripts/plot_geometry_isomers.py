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
                'DMSO - pm6': 2.973747,
                'DMSO - DFT': 0,
            },
            '1112': {
                'gas - xtb': -263.495225539274,
                'DMSO - xtb': -264.111998704079,
                'DMSO - pm6': 2.971152,
                'DMSO - DFT': 0,
            },
            '1122': {
                'gas - xtb': -263.497111563998,
                'DMSO - xtb': -264.113844235440,
                'DMSO - pm6': 2.970532,
                'DMSO - DFT': 0,
            },
            '1212': {
                'gas - xtb': -263.496107977462,
                'DMSO - xtb': -264.112553678762,
                'DMSO - pm6': 2.970205,
                'DMSO - DFT': 0,
            },
        },
        'cage2': {
            '1111': {
                'gas - xtb': -208.344281547373,
                'DMSO - xtb': -208.979647567300,
                'DMSO - pm6': 2.335850,
                'DMSO - DFT': 0,
            },
            '1112': {
                'gas - xtb': -208.348232129357,
                'DMSO - xtb': -208.983862936801,
                'DMSO - pm6': 2.332590,
                'DMSO - DFT': 0,
            },
            '1122': {
                'gas - xtb': -208.348822919629,
                'DMSO - xtb': -208.983973230052,
                'DMSO - pm6': 2.334376,
                'DMSO - DFT': 0,
            },
            '1212': {
                'gas - xtb': -208.347359487125,
                'DMSO - xtb': -208.982314660766,
                'DMSO - pm6': 2.335488,
                'DMSO - DFT': 0,
            },
        }
    }


def make_plot(filename, energies):
    fig, ax = plt.subplots()
    X_positions = {
        'all': (2, '1111'), '1-down': (4, '1112'),
        'syn': (6, '1122'), 'anti': (8, '1212')
    }
    cs = ['r', 'b']
    ms = ['o', 'X']
    solv = ['gas - xtb', 'DMSO - xtb']
    # solv = ['DMSO - xtb']
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
            if len(solv) > 1:
                label = f'{cage} - {s}'
            else:
                label = f'{cage}'
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


def make_cage_plot(filename, energies, method):
    fig, ax = plt.subplots()
    X_positions = {
        'all': (2, '1111'), '1-down': (4, '1112'),
        'syn': (6, '1122'), 'anti': (8, '1212')
    }
    cs = ['r', 'b']
    ms = ['o', 'X']
    # solv = ['gas - xtb', 'DMSO - xtb']
    solv = [f'DMSO - {method}']
    cages = ['cage1', 'cage2']
    for s in solv:
        for cage, m, c in zip(cages, ms, cs):
            X_values = []
            Y_values = []
            for X in X_positions:
                X_values.append(X_positions[X][0])
                Y_values.append(
                    energies[cage][X_positions[X][1]][s] * 2625.50
                )
            new_Y_values = [i-min(Y_values) for i in Y_values]
            print(s, cage, new_Y_values)
            if len(solv) > 1:
                label = f'{cage} - {s}'
            else:
                label = f'{cage}'
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
    ax.set_title(method, fontsize=16)
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
    Usage: plot_geometry_isomers.py
        """)
        sys.exit()

    energies = energy_data()

    filename = 'geometry_isomer_energies_xtb.pdf'
    make_plot(filename, energies)
    method = 'xtb'
    filename = f'geometry_isomer_energies_{method}_dmso.pdf'
    make_cage_plot(filename, energies, method)
    method = 'pm6'
    filename = f'geometry_isomer_energies_{method}_dmso.pdf'
    make_cage_plot(filename, energies, method)
    method = 'DFT'
    filename = f'geometry_isomer_energies_{method}_dmso.pdf'
    make_cage_plot(filename, energies, method)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
