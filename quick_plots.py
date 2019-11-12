#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build and analyse all isomers MOCs with unsymm ligands.

Author: Andrew Tarzia

Date Created: 05 Nov 2019

"""

import sys
import pandas as pd

import cage_analysis as CA


def main():
    if (not len(sys.argv) == 1):
        print(
            """
            Usage: quick_plots.py
            """
        )
        sys.exit()
    else:
        pass

    experiments = [
        'li1_lk2_li5', 'li2_lk2_li6', 'li1_lk2_li4', 'li4_lk2_li5',
    ]
    cages_cis_wins = []
    cages_not_wins = []
    lig_studied = []
    energy_preferences = []
    plane_devs = []

    data = pd.read_csv('all_cage_results.txt')

    print(data[data['plane_dev_C'] > 0.3])

    for i, row in data.iterrows():
        test1 = bool(row['stable']) is True
        test2 = bool(row['preferred']) is True
        if test1 and test2:
            cages_cis_wins.append(row['lig'])
        else:
            cages_not_wins.append(row['lig'])
        lig_studied.append(row['lig'])
        plane_devs.append(float(row['plane_dev_C']))
        energies = [
            float(row['energy_A']),
            float(row['energy_B']),
            float(row['energy_C']),
            float(row['energy_D']),
        ]
        if energies[2] != 0:
            energy_preferences.append(-energies[2])
        else:
            min_diff = min([
                i for i in energies if i != 0
            ])
            energy_preferences.append(min_diff)

    CA.plot_energetics(
        lig_studied,
        experiments,
        cages_cis_wins,
        cages_not_wins,
        energy_preferences
    )

    CA.plot_plane_devs(
        lig_studied,
        experiments,
        cages_cis_wins,
        cages_not_wins,
        plane_devs
    )


if __name__ == "__main__":
    main()
