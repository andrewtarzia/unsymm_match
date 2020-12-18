#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot parity of plane deviation and order parameter.

Author: Andrew Tarzia

Date Created: 18 Dec 2020

"""

import pandas as pd
import matplotlib.pyplot as plt

from atools import colors_i_like


def main():

    cage_data = pd.read_csv('all_cage_results.txt')

    mA = 'o'
    cA = colors_i_like()[4]
    mB = 's'
    cB = colors_i_like()[3]
    mC = 'X'
    cC = colors_i_like()[0]
    mD = 'D'
    cD = colors_i_like()[2]

    X_A = []
    Y_A = []
    X_B = []
    Y_B = []
    X_C = []
    Y_C = []
    X_D = []
    Y_D = []
    for i, row in cage_data.iterrows():
        X_A.append(float(row['sqpl_op_A']))
        Y_A.append(float(row['plane_dev_A']))
        X_B.append(float(row['sqpl_op_B']))
        Y_B.append(float(row['plane_dev_B']))
        X_C.append(float(row['sqpl_op_C']))
        Y_C.append(float(row['plane_dev_C']))
        X_D.append(float(row['sqpl_op_D']))
        Y_D.append(float(row['plane_dev_D']))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        X_A,
        Y_A,
        c=cA,
        edgecolors='k',
        marker=mA,
        alpha=1,
        s=60,
        label='A isomer'
    )
    ax.scatter(
        X_B,
        Y_B,
        c=cB,
        edgecolors='k',
        marker=mB,
        alpha=1,
        s=60,
        label='B isomer'
    )
    ax.scatter(
        X_C,
        Y_C,
        c=cC,
        edgecolors='k',
        marker=mC,
        alpha=1,
        s=60,
        label='C isomer'
    )
    ax.scatter(
        X_D,
        Y_D,
        c=cD,
        edgecolors='k',
        marker=mD,
        alpha=1,
        s=60,
        label='D isomer'
    )
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'$q_{\mathrm{sqp,min}}$', fontsize=16)
    ax.set_ylabel(
        r'max. plane deviation [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        'all_cages_geom_parity.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


if __name__ == "__main__":
    main()
