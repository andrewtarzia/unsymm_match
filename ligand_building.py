#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for building ligands.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

from os.path import exists
import stk


def ligands():
    ligand_smiles = {
        'li1': 'BrC#Cc1cccc2cnccc12',
        'li2': 'BrC#Cc1cccc2ccncc12',
        'li3': 'BrC#Cc1cccc2ncccc12',
        'li4': 'BrC#Cc1cccnc1',
        'li5': 'Brc1cccnc1',
        'li6': 'Brc1ccncc1',
        # 'li7': 'Brc1cccc2cnccc12',
        # 'li8': 'Brc1cccc2ccncc12',
        # 'li9': 'Brc1cccc2ncccc12',
    }

    return ligand_smiles


def linkers():
    linker_smiles = {
        'lk1': 'Brc1ccc2[nH]c3ccc(Br)cc3c2c1',
        'lk2': 'Brc1cccc(Br)c1',
        'lk3': 'Brc1ccc(Br)cc1',
        'lk4': 'Brc1ccc2ccc(Br)cc2c1',
    }

    return linker_smiles


def build_linker(lig1_smiles, lig2_smiles, linker_smiles, name):
    """
    Build a linker from 3 components.

    Parameters
    ----------

    Returns
    -------

    """

    opt_file = f'{name}_opt.mol'

    # Build polymer if not already done.
    if exists(opt_file):
        molecule = stk.BuildingBlock.init_from_file(
            path=opt_file,
            functional_groups=['pyridine_N_metal']
        )
    else:
        print(f'building {name}')
        # Build BuildingBlocks.
        lig1 = stk.BuildingBlock(lig1_smiles, ['bromine'])
        lig2 = stk.BuildingBlock(lig2_smiles, ['bromine'])
        link = stk.BuildingBlock(linker_smiles, ['bromine'])

        # Build polymer.
        molecule = stk.ConstructedMolecule(
            building_blocks=[lig1, link, lig2],
            topology_graph=stk.polymer.Linear(
                repeating_unit='ABC',
                num_repeating_units=1,
                orientations=(0, 0, 0),
                num_processes=1
            )
        )

        # Optimise with ETKDG.
        etkdg = stk.ETKDG(use_cache=False)
        etkdg.optimize(molecule)

        # Initialise as building block.
        molecule = stk.BuildingBlock.init_from_molecule(
            mol=molecule,
            functional_groups=['pyridine_N_metal']
        )
        # Save file.
        molecule.write(opt_file)

    return molecule
