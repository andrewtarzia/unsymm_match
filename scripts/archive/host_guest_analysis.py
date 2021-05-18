#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run simple host-guest anaylsis of selected cages.

Author: Andrew Tarzia

Date Created: 10 Sep 2020

"""

import sys
from os.path import exists, join
from itertools import product

import stk
import stko

from ..utilities import MOC_xtb_opt


def build_guests():

    guest_smiles = {
        'oF': 'C1=CC=C(C(=C1)C=O)F',
        'oB': 'C1=CC=C(C(=C1)C=O)Br',
        'oC': 'C1=CC=C(C(=C1)C=O)Cl',
        'pF': 'C1=CC(=CC=C1C=O)F',
        'pB': 'C1=CC(=CC=C1C=O)Br',
        'pC': 'C1=CC(=CC=C1C=O)Cl',
        'mF': 'C1=CC(=CC(=C1)F)C=O',
        'mB': 'C1=CC(=CC(=C1)Br)C=O',
        'mC': 'C1=CC(=CC(=C1)Cl)C=O',
    }

    guests = {}
    for guest in guest_smiles:
        opt_file = f'{guest}_unopt.mol'

        if exists(opt_file):
            guests[guest] = stk.BuildingBlock.init_from_file(opt_file)
        else:
            print(f'>> optimising {guest}')
            # Build BuildingBlocks.
            molecule = stk.BuildingBlock(guest_smiles[guest])

            # Optimise with ETKDG.
            etkdg = stko.ETKDG()
            molecule = etkdg.optimize(molecule)

            # Save file.
            molecule.write(opt_file)

            # Initialise as building block.
            guests[guest] = stk.BuildingBlock.init_from_molecule(
                molecule=molecule,
            )

    return guests


def build_complexes(hosts, guests):

    complexes = {}
    for host_name, guest_name in product(hosts, guests):
        print(f'>> building complex of {host_name} and {guest_name}')
        name = f'{host_name}_g{guest_name}'
        out_file = f'{name}_unopt.mol'

        host = hosts[host_name]
        guest = guests[guest_name]

        # pd_pd_vector = [0, 0, 0]
        # fg_fg_vector = [0, 0, 0]

        complex = stk.ConstructedMolecule(
            topology_graph=stk.host_guest.Complex(
                host=host,
                guest=guest,
                # guest_start=guest.get_direction(),
                # guest_target=pd_pd_vector,
            ),
        )
        # Save file.
        complex.write(out_file)

        # Initialise as building block.
        complexes[name] = complex

    return complexes


def load_best_conformer(molecules, run_directory):

    opt_molecules = {}
    for molecule in molecules:
        print(f'>> collecting best conformer of {molecule}')
        unop_name = f'{molecule}_unopt'
        op_name = f'{molecule}_best'
        expected_out_file = (
            f"{join(run_directory, join(unop_name, 'crest_best.xyz'))}"
        )
        if not exists(expected_out_file):
            input(f'>>>> warning - missing {expected_out_file}')
            continue

        opt_molecules[molecule] = (
            molecules[molecule].with_structure_from_file(
                expected_out_file
            )
        )
        opt_molecules[molecule].write(f'{op_name}.mol')

    return opt_molecules


def opt_best_conformer(molecules, charges):

    opt_molecules = {}
    for molecule in molecules:
        op_name = f'{molecule}_opt'
        mol = molecules[molecule]
        if exists(op_name):
            mol = mol.with_structure_from_file(op_name)
        else:
            print(f'>> optimising best conformer of {molecule}')
            mol = MOC_xtb_opt(
                mol,
                molecule,
                nc=1,
                free_e=0,
                charge=charges[molecule],
                opt_level='extreme',
                etemp=300,
                solvent=('dmso', 'verytight')
            )

            mol.write(f'{op_name}.mol')

        opt_molecules[molecule] = mol

    return opt_molecules


def calculate_crest_conformer_energies(molecules, charges, suffix):

    energies = {}
    for molecule in molecules:
        energy_file = f'{molecule}_{suffix}.ey'
        if exists(energy_file):
            # Read .ey file.
            with open(energy_file, 'r') as f:
                lines = f.readlines()
            for line in lines:
                energy = float(line.rstrip())
                break
        else:
            print(f'>> calculating energy of {molecule}')
            xtb_energy = stko.XTBEnergy(
                xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
                output_dir=f'ey_{molecule}',
                num_cores=1,
                charge=charges[molecule],
                num_unpaired_electrons=0,
                electronic_temperature=300,
                unlimited_memory=True,
                solvent='DMSO',
                solvent_grid='verytight'
            )
            energy = xtb_energy.get_energy(molecules[molecule])
            # Save to .ey file.
            with open(energy_file, 'w') as f:
                f.write(f'{energy}\n')

        # to kJ/mol.
        energies[molecule] = energy * 2625.5

    return energies


def main():
    if (not len(sys.argv) == 3):
        print(
            """
            Usage: host_guest_analysis.py cage_directory
                cage_directory (str) - directory with optimised cages.

                crest_run_dir (str) - directory where crest jobs were
                manually run.

            """
        )
        sys.exit()
    else:
        cage_directory = sys.argv[1]
        run_directory = sys.argv[2]

    selected_cages = ['5A3', '5A1']
    cages = {
        i: stk.BuildingBlock.init_from_file(
            join(cage_directory, f'{i}_C_optc.mol')
        )
        for i in selected_cages
    }

    guests = build_guests()
    guest_charges = {i: 0 for i in guests}

    complex_structures = build_complexes(hosts=cages, guests=guests)
    cage_charges = {i: 4 for i in cages}
    complex_charges = {i: 4 for i in complex_structures}

    # Collect best conformers.
    best_guest_structures = load_best_conformer(
        molecules=guests,
        run_directory=run_directory,
    )
    best_complex_structures = load_best_conformer(
        molecules=complex_structures,
        run_directory=run_directory,
    )

    # Further optimise best conformers.
    opt_guest_structures = opt_best_conformer(
        molecules=best_guest_structures,
        charges=guest_charges,
    )
    opt_complex_structures = opt_best_conformer(
        molecules=best_complex_structures,
        charges=complex_charges,
    )

    # Calculate xtb conformer energies.
    best_host_energies = calculate_crest_conformer_energies(
        molecules=cages,
        charges=cage_charges,
        suffix='best',
    )
    best_complex_energies = calculate_crest_conformer_energies(
        molecules=best_complex_structures,
        charges=complex_charges,
        suffix='best',
    )
    best_guest_energies = calculate_crest_conformer_energies(
        molecules=best_guest_structures,
        charges=guest_charges,
        suffix='best',
    )
    opt_host_energies = calculate_crest_conformer_energies(
        molecules=cages,
        charges=cage_charges,
        suffix='opt',
    )
    opt_complex_energies = calculate_crest_conformer_energies(
        molecules=opt_complex_structures,
        charges=complex_charges,
        suffix='opt',
    )
    opt_guest_energies = calculate_crest_conformer_energies(
        molecules=opt_guest_structures,
        charges=guest_charges,
        suffix='opt',
    )

    print(best_host_energies, opt_host_energies)
    print(best_complex_energies, opt_complex_energies)
    print(best_guest_energies, opt_guest_energies)

    best_binding_energies = {}
    for complex in best_complex_energies:
        host, guest = complex.split('_g')
        host_ey = best_host_energies[host]
        guest_ey = best_guest_energies[guest]
        best_binding_energies[complex] = (
            best_complex_energies[complex] - (host_ey + guest_ey)
        )

    print(best_binding_energies)

    opt_binding_energies = {}
    for complex in opt_complex_energies:
        host, guest = complex.split('_g')
        host_ey = opt_host_energies[host]
        guest_ey = opt_guest_energies[guest]
        opt_binding_energies[complex] = (
            opt_complex_energies[complex] - (host_ey + guest_ey)
        )

    print(opt_binding_energies)

    sys.exit()


if __name__ == "__main__":
    main()
