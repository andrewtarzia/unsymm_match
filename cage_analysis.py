#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for analysing cages.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

from os.path import exists
import numpy as np

import stko
import atools

from plotting import isomer_plot


def get_energy(name, cage, solvent=None):
    """
    Get xTB energy of a cage.

    Parameters
    ----------

    Returns
    -------

    """

    energy_file = f'{name}_optc.ey'
    if exists(energy_file):
        # Read .ey file.
        with open(energy_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            energy = float(line.rstrip())
            break
    else:
        print(f'calculating energy of {name}')
        # Extract energy.
        if solvent is None:
            solvent_str = None
            solvent_grid = 'normal'
        else:
            solvent_str, solvent_grid = solvent
        xtb_energy = stko.XTBEnergy(
            xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
            output_dir=f'calculations_output/cage_ey_{name}',
            num_cores=6,
            charge=4,
            num_unpaired_electrons=0,
            electronic_temperature=300,
            unlimited_memory=True,
            solvent=solvent_str,
            solvent_grid=solvent_grid
        )
        energy = xtb_energy.get_energy(cage)
        # Save to .ey file.
        with open(energy_file, 'w') as f:
            f.write(f'{energy}\n')

    return energy*2625.5


def get_cage_energies(name, cages):
    """
    Get xTB energies of all four isomers.

    Parameters
    ----------

    Returns
    -------

    """

    def horiz_lines():
        # These lines are the energy difference between the cis and
        # next most stable isomer from xtb after RDKIT opt.
        # expts = {
        #     '5D1': (17, 'k'),
        #     '4D2': (27, 'k'),
        #     '5D3': (8, 'k'),
        #     '3D1': (1.4, 'r')
        # }
        # lines = [
        #     expts['5D1'][0],
        #     expts['4D2'][0],
        #     expts['5D3'][0],
        #     expts['3D1'][0],
        # ]
        # styles = [
        #     expts['5D1'][1],
        #     expts['4D2'][1],
        #     expts['5D3'][1],
        #     expts['3D1'][1],
        # ]
        higlights = {
            'threshold': (6, 'k', 3),
            'zero': (0, 'k', 3),
        }
        lines = [
            higlights['threshold'][0],
            higlights['zero'][0],
        ]
        styles = [
            (higlights['threshold'][1], higlights['threshold'][2]),
            (higlights['zero'][1], higlights['zero'][2]),
        ]
        return lines, styles

    energies = {'A': None, 'B': None, 'C': None, 'D': None}

    for iso in energies:
        name_ = f'{name}_{iso}'
        energies[iso] = get_energy(
            name=name_,
            cage=cages[iso],
            solvent=('DMSO', 'verytight')
        )

    min_energy = min(energies.values())
    energies = {
        i: energies[i]-min_energy for i in energies
    }

    isomer_plot(
        dictionary=energies,
        file_name=f'{name}_energies_plot.pdf',
        ytitle=r'relative energies [kJmol$^{-1}$]',
        ylim=(-5, 50),
        horiz=horiz_lines()
    )
    return energies


def calculate_ligand_distortion(
    mol,
    cage_name,
    free_ligand_name,
    metal_atom_nos=None
):
    """
    Calculate ligand distorion of ligands in mol.

    Strain energy definition:
        strain energy =
            E(ligand extracted from cage) -
            E(lowest energy conformer of ligand)

    Parameters
    ----------
    mol : :class:`stk.Molecule`
        Molecules whose ligands you want to analyse.

    cage_name : :class:`str`
        Name of cage associated with ligand.

    free_ligand_name : :class:`str`
        Name of ligand.

    metal_atom_nos : :class:`iterable` of :class:`int`
        The atomic number of metal atoms to remove from structure.

    Returns
    -------
    distortions : :class:`tuple` of :class:`dicts`
        Dictionaries containing the ligand distortions of each ligand
        found in the cage.
            Distortions calculated depends on input:
                Average change in NN distance from lowest energy
                conformer.
                Average change in bite angle from lowest energy
                conformer.
                Strain energy from lowest energy conformer.

    """

    org_ligs, smiles_keys = atools.get_organic_linkers(
        cage=mol,
        metal_atom_nos=metal_atom_nos,
        file_prefix=f'{cage_name}_sg'
    )

    atools.get_lowest_energy_conformers(
        org_ligs=org_ligs,
        smiles_keys=smiles_keys,
        file_prefix=f'{free_ligand_name}_sg'
    )

    deltann_dist_dict = atools.calculate_deltann_distance(
        org_ligs=org_ligs,
        smiles_keys=smiles_keys,
        fg_factory=[
            atools.AromaticCNCFactory(),
            atools.AromaticCNNFactory()
        ],
        file_prefix=f'{free_ligand_name}_sg'
    )
    deltaangle_dist_dict = atools.calculate_deltaangle_distance(
        org_ligs=org_ligs,
        smiles_keys=smiles_keys,
        fg_factory=[
            atools.AromaticCNCFactory(),
            atools.AromaticCNNFactory()
        ],
        file_prefix=f'{free_ligand_name}_sg'
    )
    lse_dict = atools.calculate_ligand_SE(
        org_ligs=org_ligs,
        smiles_keys=smiles_keys,
        output_json=f'{cage_name}_lse.json',
        file_prefix=f'{free_ligand_name}_sg'
    )

    NN_avg_cage_min_free = np.average(list(
        deltann_dist_dict.values()
    ))
    bite_avg_cage_min_free = np.average(list(
        deltaangle_dist_dict.values()
    ))
    sum_strain_energy = sum(list(lse_dict.values()))

    distortions = (
        NN_avg_cage_min_free,
        bite_avg_cage_min_free,
        sum_strain_energy
    )

    return distortions


def get_ligand_distortion(name, cages, NN_dists, bites_dist):
    """
    Get ligand distortions compared to free ligand.

    Parameters
    ----------

    Returns
    -------

    """

    l_distortions = {
        'bite_angle': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'deltabites',
            r'average($\Delta$bite angle) [degrees]',
            (0, 100)
        ),
        'NN_dist': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'deltaNNs',
            r'average($\Delta$NN distance) [$\mathrm{\AA}$]',
            (0, 10)
        ),
        'sum_strain': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'sum_strain',
            r'sum(ligand strain energy) [kJmol$^{-1}$]',
            (-10, 300)
        ),
    }
    for iso in cages:
        name_ = f'{name}_{iso}'
        cage = cages[iso]
        results = calculate_ligand_distortion(
            mol=cage,
            cage_name=name_,
            free_ligand_name=name,
            metal_atom_nos=(46, )
        )
        NN_change, bite_change, sum_strain = results
        l_distortions['NN_dist'][0][iso] = NN_change
        l_distortions['bite_angle'][0][iso] = bite_change
        l_distortions['sum_strain'][0][iso] = sum_strain

    for i in l_distortions:
        isomer_plot(
            dictionary=l_distortions[i][0],
            file_name=f'{name}_{l_distortions[i][1]}_plot.pdf',
            ytitle=l_distortions[i][2],
            ylim=l_distortions[i][3]
        )
    return l_distortions


def get_metal_centre_distortion(name, cages):
    """
    Get metal centre distortions.

    Parameters
    ----------

    Returns
    -------

    """

    m_distortions = {
        'bond_lengths': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'bonds',
            r'mean bond length [$\mathrm{\AA}$]',
            (1, 3)
        ),
        'angles': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'angles',
            'mean N-Pd-N angle [degrees]',
            (0, 100)
        ),
        'torsions': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'torsions',
            'mean N-N-N-N torsion [degrees]',
            (0, 100)
        ),
        'plane_dev': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'planes',
            r'max. sum of plane deviation [$\mathrm{\AA}$]',
            (0, 2)
        ),
        'min_q4_op': (
            {'A': None, 'B': None, 'C': None, 'D': None},
            'q4op',
            r'min. $q_{\mathrm{sp}}$',
            (0, 1)
        )
    }

    for iso in cages:
        cage = cages[iso]
        results = atools.get_square_planar_distortion(
            mol=cage,
            metal=46,
            bonder=7
        )
        order_results = atools.get_order_values(
            mol=cage,
            metal=46
        )

        for measure in results:
            if measure == 'plane_dev':
                m_distortions[measure][0][iso] = max(results[measure])
            elif measure not in m_distortions:
                continue
            else:
                m_distortions[measure][0][iso] = np.mean(
                    results[measure]
                )

        m_distortions['min_q4_op'][0][iso] = (
            order_results['sq_plan']['min']
        )

    for i in m_distortions:
        isomer_plot(
            dictionary=m_distortions[i][0],
            file_name=f'{name}_{m_distortions[i][1]}_plot.pdf',
            ytitle=m_distortions[i][2],
            ylim=m_distortions[i][3]
        )
    return m_distortions


def check_stability(l_distortions, m_distortions):
    """
    Check if cis isomer is stable based on distortions.

    Parameters
    ----------

    Returns
    -------

    """

    def fail_NN(value):
        return value > 1

    def fail_plane_dev(value):
        return value > 0.3

    def fail_op(value):
        return value < 0.95

    def fail_sum_strain(value):
        return value < 0.7

    checks = {
        'bond_lengths': None,
        'angles': None,
        'torsions': None,
        'plane_dev': fail_plane_dev,
        'bite_angle': None,
        'NN_dist': None,
        'min_q4_op': fail_op,
        'sum_strain': None,
    }

    for i in l_distortions:
        # No test for this measure.
        if checks[i] is None:
            continue
        check = checks[i](l_distortions[i][0]['C'])
        print(
            f"> {i}: {round(l_distortions[i][0]['C'], 4)}."
            f" C isomer unstable? {check}"
        )
        if checks[i](l_distortions[i][0]['C']):
            return False

    for i in m_distortions:
        # No test for this measure.
        if checks[i] is None:
            continue
        check = checks[i](m_distortions[i][0]['C'])
        print(
            f"> {i}: {round(m_distortions[i][0]['C'], 4)}."
            f" C isomer unstable? {check}"
        )
        if check:
            return False

    return True


def check_preference(energies, energy_cutoff):
    """
    Check if cis isomer is preferred based on relative energetics.

    Parameters
    ----------

    Returns
    -------

    """

    print(f'>>> relative isomer energies {energies}')

    if energies['C'] == 0:
        energy_sep = min([
            energies[i] for i in energies if energies[i] != 0
        ])
        print(
            f"> energy test: C: {round(energies['C'], 4)}, "
            f'next: {round(energy_sep, 4)} kJmol$^{-1}$. Less than '
            f'threshold? {energy_sep < energy_cutoff}.'
        )
        if energy_sep < energy_cutoff:
            return False, energy_sep
    else:
        print(
            f"> energy test failed: C is {round(energies['C'], 4)} "
            'kJmol$^{-1}$ less stable'
        )
        return False, -energies['C']

    return True, energy_sep
