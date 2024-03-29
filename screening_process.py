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
from pandas import read_csv

import ligand_building as LB
import ligand_analysis as LA
import cage_building as CB
import cage_analysis as CA
import plotting as PL
import utilities


def build_all_ligands():
    """
    Build all ligands.

    Returns
    -------
    ligands : :class:`dict` of :class:`stk.BuildingBlock`
        Keys are names of ligands, value contains ligand
        building block.

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
            name = f"{lig1}{link}{lig2}"
            ligand = LB.build_linker(
                lig1_smiles=lig1_smiles,
                lig2_smiles=lig2_smiles,
                linker_smiles=link_smiles,
                name=name
            )
            ligands[name] = ligand

    print(f'{count} ligands built -> {count*4} cages to build.')
    return ligands


def analyse_all_ligands(params, ligands):
    """
    Analyse all ligands.

    Parameters
    ----------
    params: :class:`dict`
        Parameters of analysis - sets number of conformers to check
        geometries of.

    ligands : :class:`dict` of :class:`stk.BuildingBlock`
        Keys are names of ligands, value contains ligand
        building block.

    Returns
    -------
    ligands : :class:`dict` of :class:`tuple`
        Keys are names of ligands, tuple in value contains ligand
        building block.

    """

    for ligand in ligands:
        # if ligand not in experiments:
        #     continue
        name = ligand
        molecule = ligands[ligand]
        if exists(f'{name}_opt_chosen.mol'):
            with open(name+'_bites_dists.txt', 'r') as f:
                lines = f.readlines()
            bites_dist = [float(i.rstrip()) for i in lines]
            with open(name+'_NNs_dists.txt', 'r') as f:
                lines = f.readlines()
            NN_dist = [float(i.rstrip()) for i in lines]
            molecule = molecule.with_structure_from_file(
                f'{name}_opt_chosen.mol'
            )

        else:
            print(f'doing {name}')
            # Get conformers.
            confs, cids = LA.get_conformers(
                molecule,
                N=int(params['N'])
            )

            # Calculate and save their flexibility.
            bites_dist = LA.calc_bite_flexibility(
                molecule,
                confs,
                cids,
                name
            )

            NN_dist = LA.calc_NN_flexibility(
                molecule,
                confs,
                cids,
                name
            )

            # Plot their flexibility.
            bites_dist = LA.plot_bite_flexibility(bites_dist, name)
            NN_dist = LA.plot_NN_flexibility(NN_dist, name)

            # Select conformer with binding groups pointing the right
            # way.
            # Update molecule.
            molecule = LA.select_conformer(molecule, confs, cids, name)
            molecule.write(f'{name}_opt_chosen.mol')
            molecule.write(f'{name}_opt_chosen.xyz')
        ligands[ligand] = (molecule, bites_dist, NN_dist)

    return ligands


def build_all_cages(ligands):
    """
    Build all cages.

    Parameters
    ----------
    ligands : :class:`dict` of :class:`tuple`
        Keys are names of ligands, tuple in value contains ligand
        building block.

    Returns
    -------
    :class:`dict`
        Dictionary of cage isomers, where keys are the name of the cage
        set.

    """

    all_cage_sets = {}
    for ligand in ligands:
        # if ligand not in experiments:
        #     continue
        name = ligand
        molecule = ligands[ligand][0]
        all_cage_sets[name] = CB.build_cage_isomers(
            name=name,
            ligand=molecule,
        )

    return all_cage_sets


def analyse_all_cages(all_cage_sets, read_data):
    """
    Analyse all cages, produces many plots.

    Parameters
    ----------
    all_cage_sets : :class:`dict`
        Dictionary of cage isomers, where keys are the name of the cage
        set.

    read_data : :class:`bool`
        `True` to read already calculated data, `False` to recalculate.

    """

    cages_cis_wins = []
    cages_not_wins = []
    lig_studied = []
    energy_preferences = []
    plane_devs = []
    sqpl_ops = []
    lses = []

    experiments = ['5D1', '4D2', '5D3', '3D1']

    if read_data and exists('all_cage_results.txt'):
        data = read_csv('all_cage_results.txt')

        for i, row in data.iterrows():
            lig_name = row['lig']
            stable = row['stable']
            preferred = row['preferred']
            cis_preferred_and_stable = all([stable, preferred])
            if cis_preferred_and_stable:
                cages_cis_wins.append(lig_name)
            else:
                cages_not_wins.append(lig_name)

            lig_studied.append(lig_name)
            energies = [
                float(row['energy_A']), float(row['energy_B']),
                float(row['energy_C']), float(row['energy_D']),
            ]
            if energies[2] == 0:
                energy_preferences.append(
                    min([
                        i for i in energies if i != 0
                    ])
                )
            else:
                energy_preferences.append(-energies[2])
            plane_devs.append(float(row['plane_dev_C']))
            sqpl_ops.append(float(row['sqpl_op_C']))
            lses.append(float(row['lse_C']))

    else:
        with open('all_cage_results.txt', 'w') as f:
            f.write(
                'lig,stable,preferred,'
                'plane_dev_A,plane_dev_B,plane_dev_C,plane_dev_D,'
                'sqpl_op_A,sqpl_op_B,sqpl_op_C,sqpl_op_D,'
                'lse_A,lse_B,lse_C,lse_D,'
                'energy_A,energy_B,energy_C,energy_D\n'
            )
            for lig_name in all_cage_sets:
                # if lig_name not in experiments:
                #     continue
                print('ligand:', lig_name)
                lig_studied.append(lig_name)
                cages = all_cage_sets[lig_name]
                energies = CA.get_cage_energies(lig_name, cages)
                m_distortions = CA.get_metal_centre_distortion(
                    name=lig_name,
                    cages=cages
                )
                l_distortions = CA.get_ligand_distortion(
                    name=lig_name,
                    cages=cages,
                    # bites_dist=ligands[lig_name][1],
                    # NN_dists=ligands[lig_name][2]
                    bites_dist=None,
                    NN_dists=None
                )

                # print('bites', np.mean(ligands[lig_name][1]))
                # print('bites', np.std(ligands[lig_name][1]))
                # print('nn', np.mean(ligands[lig_name][2]))
                # print('nn', np.std(ligands[lig_name][2]))
                # for i in m_distortions['plane_dev'][0]:
                #     print(i, m_distortions['plane_dev'][0][i])
                #
                # for i in l_distortions['bite_angle'][0]:
                #     print(i, l_distortions['bite_angle'][0][i])
                #
                # for i in l_distortions['NN_dist'][0]:
                #     print(i, l_distortions['NN_dist'][0][i])
                # print('---')

                stable = CA.check_stability(
                    l_distortions=l_distortions,
                    m_distortions=m_distortions
                )
                preferred, energy_sep = CA.check_preference(
                    energies,
                    energy_cutoff=6.0
                )
                energy_preferences.append(energy_sep)
                plane_devs.append(m_distortions['plane_dev'][0]['C'])
                sqpl_ops.append(m_distortions['min_q4_op'][0]['C'])
                lses.append(l_distortions['sum_strain'][0]['C'])

                cis_preferred_and_stable = all([stable, preferred])
                if cis_preferred_and_stable:
                    cages_cis_wins.append(lig_name)
                else:
                    cages_not_wins.append(lig_name)

                f.write(
                    f'{lig_name},{stable},{preferred},'
                    f"{m_distortions['plane_dev'][0]['A']},"
                    f"{m_distortions['plane_dev'][0]['B']},"
                    f"{m_distortions['plane_dev'][0]['C']},"
                    f"{m_distortions['plane_dev'][0]['D']},"
                    f"{m_distortions['min_q4_op'][0]['A']},"
                    f"{m_distortions['min_q4_op'][0]['B']},"
                    f"{m_distortions['min_q4_op'][0]['C']},"
                    f"{m_distortions['min_q4_op'][0]['D']},"
                    f"{l_distortions['sum_strain'][0]['A']},"
                    f"{l_distortions['sum_strain'][0]['B']},"
                    f"{l_distortions['sum_strain'][0]['C']},"
                    f"{l_distortions['sum_strain'][0]['D']},"
                    f"{energies['A']},{energies['B']},"
                    f"{energies['C']},{energies['D']}\n"
                )
                print('-----------------------')

    print('-----------------------------------------')
    # Plot distribution of all cages.
    total_cages = len(cages_cis_wins) + len(cages_not_wins)
    print(
        f'{len(cages_cis_wins)} cages with cis preferred '
        f'and stable of {total_cages}.'
    )
    print('-----------------------------------------')
    print('candidate cages:')
    for i in sorted(cages_cis_wins):
        print(i)

    PL.plot_isomer_distributions()

    PL.plot_all_cages_bars(
        ligands=lig_studied,
        experiments=experiments,
        cages_cis_wins=cages_cis_wins,
        cages_not_wins=cages_not_wins,
        y_value=energy_preferences,
        y_title='stability of C isomer [kJmol$^{-1}$]',
        y_bar=6.0,
        suffix='energypreference',
    )

    PL.plot_all_cages_bars(
        ligands=lig_studied,
        experiments=experiments,
        cages_cis_wins=cages_cis_wins,
        cages_not_wins=cages_not_wins,
        y_value=lses,
        y_title=r'sum(ligand strain energy) [kJmol$^{-1}$]',
        y_bar=0.0,
        suffix='lses',
    )

    PL.plot_all_cages_bars(
        ligands=lig_studied,
        experiments=experiments,
        cages_cis_wins=cages_cis_wins,
        cages_not_wins=cages_not_wins,
        y_value=plane_devs,
        y_title=r'max. plane deviation [$\mathrm{\AA}$]',
        y_bar=0.3,
        suffix='planedevs',
    )

    PL.plot_all_cages_bars(
        ligands=lig_studied,
        experiments=experiments,
        cages_cis_wins=cages_cis_wins,
        cages_not_wins=cages_not_wins,
        y_value=sqpl_ops,
        y_title=r'$q_{\mathrm{sqp,min}}$',
        y_bar=0.95,
        suffix='sqplops',
    )

    PL.plot_energetics_and_geom(
        ligands=lig_studied,
        experiments=experiments,
        cages_cis_wins=cages_cis_wins,
        cages_not_wins=cages_not_wins,
        energy_preferences=energy_preferences,
        plane_devs=plane_devs,
        sqpl_ops=sqpl_ops,
    )
    PL.plot_energetics_and_geom_3D(
        ligands=lig_studied,
        experiments=experiments,
        energy_preferences=energy_preferences,
        plane_devs=plane_devs,
        sqpl_ops=sqpl_ops,
    )

    utilities.draw_molecules(
        ligands=lig_studied,
        energy_preferences=energy_preferences,
        plane_devs=plane_devs,
        sqpl_ops=sqpl_ops,
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
    with open(file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        key, val = line.rstrip().split(':')
        params[key] = val

    return params


def main():
    if (not len(sys.argv) == 3):
        print(
            """
            Usage: screening_process.py param_file
                param_file (str) - text file with parameters.

                read_data (str) - 't' to read all_cage_results and
                do plots.

            """
        )
        sys.exit()
    else:
        params = read_params(sys.argv[1])
        read_data = True if sys.argv[2] == 't' else False

    # Build all ligands.
    ligands = build_all_ligands()

    # Analyse all ligands.
    ligands = analyse_all_ligands(params, ligands)

    # Build all cages.
    all_cage_sets = build_all_cages(ligands)

    # Analyse all cages.
    analyse_all_cages(all_cage_sets, read_data)


if __name__ == "__main__":
    main()
