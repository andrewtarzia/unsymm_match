#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build and analyse all isomers MOCs with unsymm ligands.

Author: Andrew Tarzia

Date Created: 05 Nov 2019

"""

import sys
from os.path import exists
from itertools import combinations

import ligand_building as LB
import ligand_analysis as LA
import cage_building as CB
import cage_analysis as CA


def build_all_ligands(params):
    """
    Build all cages.

    Parameters
    ----------

    Returns
    -------

    """
    ligand_smiles = LB.ligands()
    linker_smiles = LB.linkers()

    ligands = {}

    # Iterate through linkers.
    count = 0
    for link in linker_smiles:
        link_smiles = linker_smiles[link]
        for lig1, lig2 in combinations(ligand_smiles, r=2):
            count += 1
            lig1_smiles = ligand_smiles[lig1]
            lig2_smiles = ligand_smiles[lig2]
            name = f"{lig1}_{link}_{lig2}"
            ligand = LB.build_linker(
                lig1_smiles=lig1_smiles,
                lig2_smiles=lig2_smiles,
                linker_smiles=link_smiles,
                name=name
            )
            ligands[name] = ligand

    print(f'{count} ligands built')
    return ligands


def analyse_all_ligands(params, ligands):
    """
    Analyse all ligands.

    Parameters
    ----------

    Returns
    -------

    """
def build_all_cages(params):
    """
    Build all cages.

    Parameters
    ----------

    Returns
    -------

    """
    raise NotImplementedError()


def analyse_all_cages(params):
    """
    Analyse all cages.

    Parameters
    ----------

    Returns
    -------

    """
    raise NotImplementedError()


def read_params(file):
    """
    Read parametrs for screening.

    Returns
    -------
    params : :class:`dict`
        Dictionary of parameters.

    """
    params = {}
    with open(sys.argv[1], 'r') as f:
        lines = f.readlines()

    for line in lines:
        key, val = line.rstrip().split(':')
        params[key] = val

    return params


def main():
    if (not len(sys.argv) == 2):
        print(
            """
            Usage: screening_process.py param_file
                param_file (str) - text file with parameters.
            """
        )
        sys.exit()
    else:
        params = read_params(sys.argv[1])

    print(params)

    # Build all ligands.
    ligands = build_all_ligands(params)
    # Analyse all ligands.
    ligands = analyse_all_ligands(params, ligands)
    # Build all cages.
    build_all_cages(params)
    # Analyse all cages.
    analyse_all_cages(params)


if __name__ == "__main__":
    main()
