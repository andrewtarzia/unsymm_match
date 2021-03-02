#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities for scripts.

Author: Andrew Tarzia

Date Created: 02 Mar 2021

"""


import matplotlib.pyplot as plt


def bar_chart(energies, filename, x_labels, y_title, facecolor):

    xs = [1, 2, 3, 4]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(
        xs,
        energies,
        width=1,
        color=facecolor,
        edgecolor='k',
        linewidth=2
    )
    ax.set_xticks(xs)
    ax.set_xticklabels(x_labels)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand', fontsize=16)
    ax.set_ylabel(y_title, fontsize=16)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, None)

    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
    plt.close()


def get_g16_energy(
    filename,
    front_splitter,
    back_splitter,
    back_splitter2,
):

    with open(filename, 'r') as f:
        data = f.read()

    energy = data.split(front_splitter)
    energy = energy[-1].split(back_splitter)[0]
    try:
        return float(energy)  # a.u.
    except ValueError:
        if back_splitter2 is None:
            raise ValueError()
        energy = energy.split(back_splitter2)[0]
        return float(energy)  # a.u.


def collate_energies(file_list, functional):

    collated_energies = {}
    for file in file_list:
        splits = file.split('/')
        mol_name = splits[1]
        splits = mol_name.split('_')
        mol_name = splits[1].lower()+'_'+splits[2][0].upper()
        energy = get_g16_energy(
            filename=file,
            front_splitter=f'SCF Done:  E({functional}) =',
            back_splitter='A.U. after',
            back_splitter2=None,
        )
        collated_energies[mol_name] = energy

    return collated_energies
