#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot comparison of solvated and not DFT energies.

Author: Andrew Tarzia

Date Created: 24 Feb 2021

"""

import json
import matplotlib.pyplot as plt


def plot_parity(w_solv, wo_solv):

    fig, ax = plt.subplots(figsize=(8, 5))
    w_energies = []
    wo_energies = []
    for i in w_solv:
        w_energies.append(w_solv[i])
        wo_energies.append(wo_solv[i])

    ax.scatter(
        w_energies,
        wo_energies,
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=120
    )
    # ax.plot(
    #     np.linspace(-5000, -3000, 10), np.linspace(-5000, -3000, 10),
    #     c='k',
    #     lw=2,
    # )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'with solvation [kJmol$^{-1}$]', fontsize=16)
    ax.set_ylabel(r'without solvation [kJmol$^{-1}$]', fontsize=16)
    ax.set_lo
    fig.tight_layout()
    fig.savefig('solv_vs_nosolv.pdf', dpi=720, bbox_inches='tight')

    plt.close()


def main():

    solvated_file = 'collated_solv_energies.json'
    no_solvated_file = 'collated_energies.json'

    with open(solvated_file, 'r') as f:
        solv_energies = json.load(f)
    with open(no_solvated_file, 'r') as f:
        no_solv_energies = json.load(f)

    print(solv_energies, no_solv_energies)

    plot_parity(solv_energies, no_solv_energies)


if __name__ == '__main__':
    main()
