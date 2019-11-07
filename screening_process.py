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

    for ligand in ligands:
        name = ligand
        molecule = ligands[ligand]
        if exists(f'{name}_opt_chosen.mol'):
            molecule.update_from_file(f'{name}_opt_chosen.mol')
        else:
            print(f'doing {name}')
            # Get conformers.
            confs, cids = LA.get_conformers(
                molecule,
                N=int(params['N'])
            )

            # Plot their flexibility.
            LA.plot_bite_flexibility(molecule, confs, cids, name)

            # Select conformer with binding groups pointing the right
            # way.
            # Update molecule.
            molecule = LA.select_conformer(molecule, confs, cids, name)
            molecule.write(f'{name}_opt_chosen.mol')
            molecule.write(f'{name}_opt_chosen.xyz')
        ligands[ligand] = molecule

    return ligands


def build_all_cages(params, ligands):
    """
    Build all cages.

    Parameters
    ----------

    Returns
    -------

    """

    # Build metal complex.
    complex = CB.build_metal_centre()

    all_cage_sets = {}
    for ligand in ligands:
        name = ligand
        molecule = ligands[ligand]
        all_cage_sets[name] = CB.build_cage_isomers(
            name=name,
            ligand=molecule,
            complex=complex
        )

    return all_cage_sets


def analyse_all_cages(params, all_cage_sets, ligands):
    """
    Analyse all cages.

    Parameters
    ----------

    Returns
    -------

    """

    cages_cis_wins = []
    cages_not_wins = []

    for lig_name in all_cage_sets:
        cages = all_cage_sets[lig_name]
        energies = CA.get_cage_energies(lig_name, cages)
        m_distortions = CA.get_metal_centre_distortion(lig_name, cages)
        l_distortions = CA.get_ligand_distortion(lig_name, cages)

        print(energies)
        print(l_distortions)
        print(m_distortions)

        cis_preferred_and_stable = all([
            CA.check_stability(l_distortions, m_distortions),
            CA.check_preference(energies, energy_cutoff=100)
        ])

        if cis_preferred_and_stable:
            cages_cis_wins.append(lig_name)
        else:
            cages_not_wins.append(lig_name)

    # Plot distribution of all cages.
    total_cages = len(cages_cis_wins) + len(cages_not_wins)
    print(
        f'{len(cages_cis_wins)} cages with cis preffered '
        f'and stable of {total_cages}.'
    )


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
    all_cage_sets = build_all_cages(params, ligands)
    # Analyse all cages.
    analyse_all_cages(params, all_cage_sets, ligands)


if __name__ == "__main__":
    main()
