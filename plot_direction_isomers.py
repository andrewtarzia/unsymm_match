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
    Returns energy of all isomers from all methods.

    Returns
    -------
    :class:`list`
        List of building block placements to define cage.

    """
    return {
        '1111': {
            'gas_xtb': -238.866118443394,
            'DMSO_xtb': -239.479582875780,
            'MeCN_xtb': -239.467189116402,
            'gas_pm6': 2.31646400444,
            'DMSO_pm6': 1.76830763472,
            'MeCN_pm6': 1.77275320582,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
            'gas_xtb_fe': -237.870367512211,
            'DMSO_xtb_fe': -238.485080278423,
            'MeCN_xtb_fe': -238.472109775210,
            'gas_pm6_fe': 3.273279,
            'DMSO_pm6_fe': 2.731675,
            'MeCN_pm6_fe': 2.739022,
            'gas_dft_fe': -3926.086846,
            'DMSO_dft_fe': -3926.613729,
            'MeCN_dft_fe': -3926.609022,
        },
        '1-111': {
            'gas_xtb': -238.875561746922,
            'DMSO_xtb': -239.488637355229,
            'MeCN_xtb': -239.476575082323,
            'gas_pm6': 2.31404865815,
            'DMSO_pm6': 1.76542597351,
            'MeCN_pm6': 1.76984625580,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
            'gas_xtb_fe': -237.880551665970,
            'DMSO_xtb_fe': -238.494223076428,
            'MeCN_xtb_fe': -238.481808576593,
            'gas_pm6_fe': 3.270733,
            'DMSO_pm6_fe': 2.730784,
            'MeCN_pm6_fe': 2.735974,
            'gas_dft_fe': -3926.090615,
            'DMSO_dft_fe': -3926.615733,
            'MeCN_dft_fe': -3926.610609,
        },
        '11-1-1': {
            'gas_xtb': -238.876692140118,
            'DMSO_xtb': -239.488138214914,
            'MeCN_xtb': -239.476614754619,
            'gas_pm6': 2.31379365889,
            'DMSO_pm6': 1.76432325280,
            'MeCN_pm6': 1.76874795925,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
            'gas_xtb_fe': -237.881888488910,
            'DMSO_xtb_fe': -238.494608200064,
            'MeCN_xtb_fe': -238.482336849073,
            'gas_pm6_fe': 3.270733,
            'DMSO_pm6_fe': 2.732914,
            'MeCN_pm6_fe': 2.738094,
            'gas_dft_fe': -3926.091721,
            'DMSO_dft_fe': 0,
            'MeCN_dft_fe': -3926.605723,
        },
        '1-11-1': {
            'gas_xtb': -238.882446232199,
            'DMSO_xtb': -239.495093077175,
            'MeCN_xtb': -239.483584362788,
            'gas_pm6': 2.31232351962,
            'DMSO_pm6': 1.76716283445,
            'MeCN_pm6': 1.77152132698,
            'gas_dft': 0,
            'DMSO_dft': 0,
            'MeCN_dft': 0,
            'gas_xtb_fe': -237.888559096884,
            'DMSO_xtb_fe': -238.502359674851,
            'MeCN_xtb_fe': -238.490428615220,
            'gas_pm6_fe': 3.266756,
            'DMSO_pm6_fe': 2.729810,
            'MeCN_pm6_fe': 2.735474,
            'gas_dft_fe': -3926.094601,
            'DMSO_dft_fe': 0,
            'MeCN_dft_fe': 0,
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

    solv = ['gas_xtb_fe', 'DMSO_xtb_fe', 'MeCN_xtb_fe']
    filename = 'slide_3_isomer_energies_xtb_fe.pdf'
    three_plots(solv, filename, energies)

    solv = ['gas_pm6_fe', 'DMSO_pm6_fe', 'MeCN_pm6_fe']
    filename = 'slide_3_isomer_energies_pm6_fe.pdf'
    three_plots(solv, filename, energies)

    solv = ['gas_dft_fe', 'DMSO_dft_fe', 'MeCN_dft_fe']
    filename = 'slide_3_isomer_energies_dft_fe.pdf'
    three_plots(solv, filename, energies)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
