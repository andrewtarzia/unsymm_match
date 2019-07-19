#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot energy of isomers of a MOC with nonsymmetric linkers
based on geometry + steric effects.

Author: Andrew Tarzia

Date Created: 19 Jul 2019

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
        'cage_gs': {
            '1111': {
                'gas - xtb': -220.946487181184,
                'DMSO - xtb': -221.575666592182,
            },
            '1112': {
                'gas - xtb': -220.975513115836,
                'DMSO - xtb': -221.605756203057,
            },
            '1122': {
                'gas - xtb': -220.976854627122,
                'DMSO - xtb': -221.604175939634,
            },
            '1212': {
                'gas - xtb': -220.979131209469,
                'DMSO - xtb': -221.605816156431,
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
    cages = ['cage_gs']
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
    ax.set_ylim(-5, 90)
    ax.set_xticks([X_positions[i][0] for i in X_positions])
    ax.set_xticklabels(list(X_positions.keys()))
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )


def make_cage_plot(filename, energies):
    fig, ax = plt.subplots()
    X_positions = {
        'all': (2, '1111'), '1-down': (4, '1112'),
        'syn': (6, '1122'), 'anti': (8, '1212')
    }
    cs = ['r', 'b']
    ms = ['o', 'X']
    # solv = ['gas - xtb', 'DMSO - xtb']
    solv = ['DMSO - xtb']
    cages = ['cage_gs']
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
    ax.set_ylim(-5, 90)
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
    Usage: plot_geometry_steric_isomers.py
        """)
        sys.exit()

    energies = energy_data()

    filename = 'geom_steric_isomer_energies_xtb.pdf'
    make_plot(filename, energies)
    filename = 'geom_steric_isomer_energies_xtb_dmso.pdf'
    make_cage_plot(filename, energies)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
