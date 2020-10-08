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
        opt_file = f'{guest}_opt.mol'

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
        print(host_name, guest_name)
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

def main():
    if (not len(sys.argv) == 2):
        print(
            """
            Usage: host_guest_analysis.py cage_directory
                cage_directory (str) - directory with optimised cages.

            """
        )
        sys.exit()
    else:
        cage_directory = sys.argv[1]

    selected_cages = ['5A3', '5A1']
    cages = {
        i: stk.BuildingBlock.init_from_file(
            join(cage_directory, f'{i}_C_optc.mol')
        )
        for i in selected_cages
    }

    guests = build_guests()

    complex_structures = build_complexes(hosts=cages, guests=guests)


if __name__ == "__main__":
    main()
