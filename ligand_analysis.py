#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for analysing ligands.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

import sys
from os.path import exists
import numpy as np
import stk
from matplotlib.pyplot import close
from rdkit.Chem import AllChem as rdkit

import atools


def plot_bite_flexibility(molecule, confs, cids, name):
    """
    Plot the flexibility of all conformers.

    Parameters
    ----------

    Returns
    -------

    """
    bites = []
    for cid in cids:
        # Need to define the functional groups.
        new_mol = stk.BuildingBlock.init_from_molecule(
            mol=molecule,
            functional_groups=['pyridine_N_metal']
        )

        # Update stk_mol to conformer geometry.
        new_mol = atools.update_from_rdkit_conf(
            new_mol,
            confs,
            conf_id=cid
        )

        bites.append(
            atools.calculate_bite_angle(
                bb=new_mol,
                constructed=False,
            )
        )
    fig, ax = atools.histogram_plot_N(
        Y=bites, X_range=(0, 180), width=5,
        alpha=1.0,
        color=atools.colors_i_like()[1],
        edgecolor=atools.colors_i_like()[1],
        xtitle=r'bite angles [degrees]',
        N=1
    )
    fig.tight_layout()
    fig.savefig(name+'_bites_dists.pdf', dpi=720,
                bbox_inches='tight')
    close()


def get_conformers(molecule, N):
    """
    Get N conformers with ETKDG algorithm.

    Parameters
    ----------

    Returns
    -------

    """

    confs = molecule.to_rdkit_mol()
    etkdg = rdkit.ETKDG()
    etkdg.randomSeed = 1000
    cids = rdkit.EmbedMultipleConfs(
        mol=confs,
        numConfs=N,
        params=etkdg
    )

    return confs, cids


def calculate_NN_angle(bb):
    """
    Calculate the N-COM-N angle of a ditopic building block.

    This function will not work for cages built from FGs other than
    metals + pyridine_N_metal.

    Parameters
    ----------
    bb : :class:`stk.BuildingBlock`
        stk molecule to analyse.

    Returns
    -------
    angle : :class:`float`
        Angle between two bonding vectors of molecule.

    """

    fg_counts = 0
    fg_positions = []
    for fg in bb.func_groups:
        if 'pyridine_N_metal' == fg.fg_type.name:
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            N_position = bb.get_center_of_mass(
                atom_ids=fg.get_deleter_ids()
            )
            fg_positions.append(N_position)

    if fg_counts != 2:
        sys.exit(
            f'{bb} does not have 2 pyridine_N_metal functional '
            'groups.'
        )

    # Get building block COM.
    COM_position = bb.get_center_of_mass()

    # Get vectors.
    fg_vectors = [i-COM_position for i in fg_positions]

    # Calculate the angle between the two vectors.
    angle = np.degrees(atools.angle_between(*fg_vectors))
    return angle


def select_conformer(molecule, confs, cids, name):
    """
    Select and optimize a conformer with desired directionality.

    Currently:
        Best directionality will be defined by the smallest
        N-ligand centroid-N angle.

    Parameters
    ----------

    Returns
    -------

    """

    if exists(f'{name}_opt_chosen.mol'):
        molecule.update_from_file(f'{name}_opt_chosen.mol')
    else:
        print(f'getting optimal conformer of {name}')
        min_angle = 10000
        min_cid = -10
        for cid in cids:
            # Need to define the functional groups.
            new_mol = stk.BuildingBlock.init_from_molecule(
                mol=molecule,
                functional_groups=['pyridine_N_metal']
            )

            # Update stk_mol to conformer geometry.
            new_mol = atools.update_from_rdkit_conf(
                new_mol,
                confs,
                conf_id=cid
            )

            angle = calculate_NN_angle(new_mol)

            if angle < min_angle:
                min_cid = cid
                min_angle = angle
                molecule = atools.update_from_rdkit_conf(
                    molecule,
                    confs,
                    conf_id=min_cid
                )
        xtb_opt = stk.XTB(
            xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
            output_dir=f'{name}',
            gfn_version=2,
            num_cores=6,
            opt_level='verytight',
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True
        )
        xtb_opt.optimize(mol=molecule)

    return molecule
