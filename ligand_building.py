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
import stko
from atools import AromaticCNCFactory, AromaticCNNFactory


def ligands():
    ligand_smiles = {
        'li1': 'BrC#Cc1cccc2cnccc12',
        'li2': 'BrC#Cc1cccc2ccncc12',
        'li3': 'BrC#Cc1cccc2ncccc12',
        'li4': 'BrC#Cc1cccnc1',
        'li5': 'Brc1cccnc1',
        'li6': 'Brc1ccncc1',
        # 'li7': 'Cc1cn(Br)nn1',  # Removed on 02/09/20
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
            functional_groups=[
                AromaticCNCFactory(),
                AromaticCNNFactory()
            ],
        )
    else:
        print(f'>> building {name}')
        # Build BuildingBlocks.
        lig1 = stk.BuildingBlock(lig1_smiles, [stk.BromoFactory()])
        lig2 = stk.BuildingBlock(lig2_smiles, [stk.BromoFactory()])
        link = stk.BuildingBlock(linker_smiles, [stk.BromoFactory()])

        # Build polymer.
        molecule = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(lig1, link, lig2),
                repeating_unit='ABC',
                num_repeating_units=1,
                orientations=(0, 0, 0),
                num_processes=1
            )
        )

        # Optimise with ETKDG.
        etkdg = stko.ETKDG()
        molecule = etkdg.optimize(molecule)

        # Initialise as building block.
        molecule = stk.BuildingBlock.init_from_molecule(
            molecule=molecule,
            functional_groups=[
                AromaticCNCFactory(),
                AromaticCNNFactory()
            ],
        )
        # Save file.
        molecule.write(opt_file)

    return molecule
