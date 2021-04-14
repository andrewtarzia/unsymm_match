#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run simple cage analysis for final figure of paper.

Author: Andrew Tarzia

Date Created: 05 Dec 2020

"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import exists
from os import system
import pywindow as pw
import json

import stk

from atools import colors_i_like, angle_between, get_atom_distance


def get_pore_size(name, molecule):

    if not exists(f'{name}_pw.json'):
        print(f'....analyzing porosity of {name}')
        # Load cage into pywindow.
        molecule.write('temp.xyz')
        pw_cage = pw.MolecularSystem.load_file('temp.xyz')
        pw_cage_mol = pw_cage.system_to_molecule()
        system('rm temp.xyz')

        # Calculate pore size.
        try:
            pw_cage_mol.calculate_pore_diameter_opt()
            pw_cage_mol.calculate_centre_of_mass()

        except ValueError:
            print(f'Warning: {name} failed pyWindow calculation!')
            # Handle failure.
            pw_cage_mol.properties = {
                'centre_of_mass': [0, 0, 0],
                'pore_diameter_opt': {
                    'atom_1': 0,
                    'centre_of_mass': [0, 0, 0],
                    'diameter': -1
                },
            }

        # Save files.
        molecule.write(f'{name}_pw.xyz')
        with open(f'{name}_pw.xyz', 'r') as f:
            lines = f.readlines()
        lines[0] = f'{len(list(molecule.get_atoms()))+2}\n'
        x, y, z = pw_cage_mol.properties['pore_diameter_opt'][
            'centre_of_mass'
        ]
        lines.append(
            f'Xe {round(x, 5)} {round(y, 5)} {round(z, 5)}\n'
        )
        x, y, z = next(molecule.get_atomic_positions(atom_ids=(
            pw_cage_mol.properties['pore_diameter_opt'][
                'atom_1'
            ],
        )))
        lines.append(
            f'Kr {round(x, 5)} {round(y, 5)} {round(z, 5)}\n'
        )
        with open(f'{name}_pw.xyz', 'w') as f:
            f.write(''.join(lines))
        pw_cage_mol.dump_properties_json(f'{name}_pw.json')

    # Get data.
    with open(f'{name}_pw.json', 'r') as f:
        data = json.load(f)

    return data['pore_diameter_opt']['diameter']


def get_distance(molecule):

    # Get 1 Pd atom and its neighbours.
    metal_atoms = []
    for atom in molecule.get_atoms():
        if atom.get_atomic_number() == 46:
            metal_atoms.append(atom)
        if len(metal_atoms) == 2:
            break

    return get_atom_distance(
        molecule=molecule,
        atom1_id=metal_atoms[0].get_id(),
        atom2_id=metal_atoms[1].get_id(),
    )


def get_deviation(molecule):

    # Get 1 Pd atom and its neighbours.
    metal_atoms = []
    for atom in molecule.get_atoms():
        if atom.get_atomic_number() == 46:
            metal_atoms.append(atom)
        if len(metal_atoms) == 2:
            break
    chosen_atom = metal_atoms[0]
    second_atom = metal_atoms[1]

    # Get plane defined by them.
    plane_ids = [chosen_atom.get_id()]
    for bond in molecule.get_bonds():
        a1 = bond.get_atom1()
        a2 = bond.get_atom2()
        if chosen_atom.get_id() == a1.get_id():
            plane_ids.append(a2.get_id())
        elif chosen_atom.get_id() == a2.get_id():
            plane_ids.append(a1.get_id())

    normal = molecule.get_plane_normal(atom_ids=plane_ids)
    # Get Pd-Pd vector.
    pos_mat = molecule.get_position_matrix()
    vector = (
        pos_mat[second_atom.get_id()] - pos_mat[chosen_atom.get_id()]
    )

    # Project Pd-Pd vector onto plane.
    # Apply from https://www.geeksforgeeks.org/vector-projection-
    # using-python/
    projection = vector - (
        np.dot(vector, normal)/(np.linalg.norm(normal)**2)
    )*normal
    projection_norm = np.linalg.norm(projection)

    # Get angle between plane and Pd-Pd vector.
    angle_between_vectors = np.degrees(
        angle_between(v1=vector, v2=normal)
    )

    return angle_between_vectors, projection_norm


def plot(
    x, y, xlabel, ylabel, xlim, ylim,
    expt_x, expt_y, sele_x, sele_y, filename,
    exam_x, exam_y, example_cases,
):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        exam_x,
        exam_y,
        c='r',
        edgecolors='r',
        marker='o',
        alpha=1,
        s=260,
    )
    ax.scatter(
        x,
        y,
        c='gold',
        edgecolors='k',
        marker='o',
        alpha=1,
        s=120,
        label='this work'
    )
    ax.scatter(
        sele_x,
        sele_y,
        c=colors_i_like()[4],
        edgecolors='k',
        marker='X',
        alpha=1,
        s=120,
        label='selected cage ligands'
    )
    ax.scatter(
        expt_x,
        expt_y,
        c=colors_i_like()[3],
        edgecolors='k',
        marker='P',
        alpha=1,
        s=120,
        label='published examples'
    )

    # ax.axhline(y=0, lw=2, linestyle='--', c='k')

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(filename, dpi=720, bbox_inches='tight')
    plt.close()


def main():
    if (not len(sys.argv) == 1):
        print(
            """
            Usage: plot_cage_properties_fig.py

            """
        )
        sys.exit()
    else:
        pass

    cage_data = pd.read_csv('all_cage_results.txt')
    experimental_ligands = ['5D3', '5D1', '4D2', '3D1']
    selected_ligands = ['4B3', '4B1', '5B4', '5A3', '5A1']
    example_cases = ['3B1', '5D2', '3C2', '5A2', '4A2', '6C1']

    cis_cages = {}
    expt_cages = {}
    sele_cages = {}
    exam_cages = {}
    for i, row in cage_data.iterrows():
        name = row['lig']
        if float(row['energy_C']) != 0:
            continue
        if float(row['sqpl_op_C']) < 0.9:
            continue
        molecule = stk.BuildingBlock.init_from_file(
            f'{name}_C_optc.mol'
        )
        c_data = {}
        c_data['pore_size'] = get_pore_size(name, molecule)
        c_data['angle'], c_data['deviation'] = get_deviation(molecule)
        c_data['distance'] = get_distance(molecule)
        if name in example_cases:
            exam_cages[name] = c_data
        if name in experimental_ligands:
            expt_cages[name] = c_data
        elif name in selected_ligands:
            sele_cages[name] = c_data
        else:
            cis_cages[name] = c_data
        if c_data['pore_size'] == -1:
            print('pywindow fail', name, c_data)
        elif c_data['pore_size'] < 3 and c_data['deviation'] > 5:
            print('small pore, high aniso', name, c_data)
        elif c_data['pore_size'] > 3 and c_data['deviation'] > 6:
            print('large pore, high aniso', name, c_data)

    plot(
        x=[cis_cages[i]['pore_size'] for i in cis_cages],
        y=[cis_cages[i]['angle'] for i in cis_cages],
        xlabel=r'pore diameter [$\mathrm{\AA}$]',
        ylabel=r'Pd angle deviation [$\mathrm{\cdot}$]',
        xlim=(0, None),
        ylim=(None, None),
        expt_x=[expt_cages[i]['pore_size'] for i in expt_cages],
        expt_y=[expt_cages[i]['angle'] for i in expt_cages],
        sele_x=[sele_cages[i]['pore_size'] for i in sele_cages],
        sele_y=[sele_cages[i]['angle'] for i in sele_cages],
        filename='all_cages_prop_angle_vpore.pdf',
        exam_x=[exam_cages[i]['pore_size'] for i in exam_cages],
        exam_y=[exam_cages[i]['angle'] for i in exam_cages],
        example_cases=example_cases,
    )
    plot(
        x=[cis_cages[i]['distance'] for i in cis_cages],
        y=[cis_cages[i]['angle'] for i in cis_cages],
        xlabel=r'Pd-Pd distance [$\mathrm{\AA}$]',
        ylabel=r'Pd angle deviation [$\mathrm{\cdot}$]',
        xlim=(None, None),
        ylim=(None, None),
        expt_x=[expt_cages[i]['distance'] for i in expt_cages],
        expt_y=[expt_cages[i]['angle'] for i in expt_cages],
        sele_x=[sele_cages[i]['distance'] for i in sele_cages],
        sele_y=[sele_cages[i]['angle'] for i in sele_cages],
        filename='all_cages_prop_angle_vdist.pdf',
        exam_x=[exam_cages[i]['distance'] for i in exam_cages],
        exam_y=[exam_cages[i]['angle'] for i in exam_cages],
        example_cases=example_cases,
    )
    plot(
        x=[cis_cages[i]['pore_size'] for i in cis_cages],
        y=[cis_cages[i]['deviation'] for i in cis_cages],
        xlabel=r'pore diameter [$\mathrm{\AA}$]',
        ylabel=r'$\Delta_{\mathrm{Pd}}$ [$\mathrm{\AA}$]',
        xlim=(0, 9),
        ylim=(None, 12),
        expt_x=[expt_cages[i]['pore_size'] for i in expt_cages],
        expt_y=[expt_cages[i]['deviation'] for i in expt_cages],
        sele_x=[sele_cages[i]['pore_size'] for i in sele_cages],
        sele_y=[sele_cages[i]['deviation'] for i in sele_cages],
        filename='all_cages_prop_dev_vpore.pdf',
        exam_x=[exam_cages[i]['pore_size'] for i in exam_cages],
        exam_y=[exam_cages[i]['deviation'] for i in exam_cages],
        example_cases=example_cases,
    )
    plot(
        x=[cis_cages[i]['distance'] for i in cis_cages],
        y=[cis_cages[i]['deviation'] for i in cis_cages],
        xlabel=r'Pd-Pd distance [$\mathrm{\AA}$]',
        ylabel=r'$\Delta_{\mathrm{Pd}}$ [$\mathrm{\AA}$]',
        xlim=(None, 15),
        ylim=(None, 12),
        expt_x=[expt_cages[i]['distance'] for i in expt_cages],
        expt_y=[expt_cages[i]['deviation'] for i in expt_cages],
        sele_x=[sele_cages[i]['distance'] for i in sele_cages],
        sele_y=[sele_cages[i]['deviation'] for i in sele_cages],
        filename='all_cages_prop_dev_vdist.pdf',
        exam_x=[exam_cages[i]['distance'] for i in exam_cages],
        exam_y=[exam_cages[i]['deviation'] for i in exam_cages],
        example_cases=example_cases,
    )

    print('--------')
    for i in cis_cages:
        print(i, cis_cages[i]['deviation'])
    print('--------')
    print(exam_cages)


if __name__ == "__main__":
    main()
