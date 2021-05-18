#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write input for ORCA 4.2.1 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

import json
import matplotlib.pyplot as plt

from utilities import replace


def plot_per_grid(edict):

    fig, ax = plt.subplots(figsize=(8, 5))
    for i in edict:

        d = edict[i]
        ax.plot(
            [int(key) for key in d],
            [
                float(d[key]['energy']) - float(d['7']['energy'])
                for key in d
            ],
            c='k',
            marker='o',
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('grid level', fontsize=16)
    ax.set_ylabel('rel. total energy [a.u.]', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        'gridlevelplot.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()
    fig, ax = plt.subplots(figsize=(8, 5))
    for i in edict:

        d = edict[i]
        ax.plot(
            [int(key) for key in d if int(key) > 3],
            [
                float(d[key]['energy']) - float(d['7']['energy'])
                for key in d if int(key) > 3
            ],
            c='k',
            marker='o',
        )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('grid level', fontsize=16)
    ax.set_ylabel('rel. total energy [a.u.]', fontsize=16)
    fig.tight_layout()
    fig.savefig(
        'gridlevelplot_zoomed.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def get_isomer(isomers, string):

    for i in isomers:
        if i in string:
            return i


def plot_isomers(edict, grid):

    isomers = ['A', 'B', 'C', 'D']
    seen_prefix = set()
    isomer_energies = {}
    for i in edict:
        pref = replace(
            string=i,
            substitutions={
                '_optc': '',
                '_A': '',
                '_B': '',
                '_C': '',
                '_D': '',
            }
        )
        isomer = get_isomer(isomers=isomers, string=i)
        if pref in seen_prefix:
            isomer_energies[pref][isomer] = (
                float(edict[i][str(grid)]['energy'])*2625.5
            )
        else:
            isomer_energies[pref] = {}
            isomer_energies[pref][isomer] = (
                float(edict[i][str(grid)]['energy'])*2625.5
            )
            seen_prefix.add(pref)

    print(isomer_energies)
    X_positions = {'A': 2, 'B': 4, 'C': 6, 'D': 8}
    for cage in isomer_energies:
        print(cage)
        cage_dict = isomer_energies[cage]
        print(cage_dict)
        min_energy = min([
            cage_dict[iso] for iso in cage_dict
        ])
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(
            [X_positions[iso] for iso in cage_dict],
            [cage_dict[iso]-min_energy for iso in cage_dict],
            c='k',
            marker='o',
        )

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('isomer', fontsize=16)
        ax.set_ylabel('rel. total energy [kJmol$^{-1}$]', fontsize=16)
        ax.set_xticks([X_positions[i] for i in X_positions])
        ax.set_xticklabels(list(X_positions.keys()))
        ax.set_ylim(-1, 50)
        fig.tight_layout()
        fig.savefig(
            f'{cage}_energy.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def main():

    grid_chosen = 6

    with open('collated_energies.json', 'r') as f:
        energies = json.load(f)

    plot_per_grid(energies)

    plot_isomers(energies, grid_chosen)


if __name__ == '__main__':
    main()
