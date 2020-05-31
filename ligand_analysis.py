#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for analysing ligands.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

from os.path import exists
import stk
import stko
from matplotlib.pyplot import close
from rdkit.Chem import AllChem as rdkit

from atools import (
    update_from_rdkit_conf,
    calculate_NN_distance,
    calculate_bite_angle,
    NPyridineFactory,
    NTriazoleFactory,
    calculate_N_COM_N_angle,
    histogram_plot_N,
    colors_i_like,
)


def calc_NN_flexibility(molecule, confs, cids, name):
    """
    Plot the flexibility of all conformers.

    Parameters
    ----------

    Returns
    -------

    """
    NNs = []
    for cid in cids:
        # Need to define the functional groups.
        new_mol = stk.BuildingBlock.init_from_molecule(
            molecule=molecule,
            functional_groups=[
                NPyridineFactory(), NTriazoleFactory()
            ]
        )

        # Update stk_mol to conformer geometry.
        new_mol = update_from_rdkit_conf(
            stk_mol=new_mol,
            rdk_mol=confs,
            conf_id=cid
        )

        NNs.append(calculate_NN_distance(bb=new_mol))

    with open(name+'_NNs_dists.txt', 'w') as f:
        for i in NNs:
            f.write(str(i)+'\n')

    return NNs


def calc_bite_flexibility(molecule, confs, cids, name):
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
            molecule=molecule,
            functional_groups=[
                NPyridineFactory(), NTriazoleFactory()
            ]
        )

        # Update stk_mol to conformer geometry.
        new_mol = update_from_rdkit_conf(
            stk_mol=new_mol,
            rdk_mol=confs,
            conf_id=cid
        )

        bites.append(calculate_bite_angle(bb=new_mol))

    with open(name+'_bites_dists.txt', 'w') as f:
        for i in bites:
            f.write(str(i)+'\n')

    return bites


def plot_NN_flexibility(NNs, name):
    """
    Plot the flexibility of all conformers.

    Parameters
    ----------

    Returns
    -------

    """
    fig, ax = histogram_plot_N(
        Y=NNs, X_range=(0, 30), width=2,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'NN distance [$\mathrm{\AA}$]',
        N=1
    )
    fig.tight_layout()
    fig.savefig(
        name+'_NNs_dists.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    close()


def plot_bite_flexibility(bites, name):
    """
    Plot the flexibility of all conformers.

    Parameters
    ----------

    Returns
    -------

    """

    fig, ax = histogram_plot_N(
        Y=bites, X_range=(0, 180), width=5,
        alpha=1.0,
        color=colors_i_like()[1],
        edgecolor=colors_i_like()[1],
        xtitle=r'bite angles [degrees]',
        N=1
    )
    fig.tight_layout()
    fig.savefig(
        name+'_bites_dists.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    close()

    return bites


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
        molecule = molecule.with_structure_from_file(
            f'{name}_opt_chosen.mol'
        )
    else:
        print(f'getting optimal conformer of {name}')
        min_angle = 10000
        min_cid = -10
        for cid in cids:
            # Need to define the functional groups.
            new_mol = stk.BuildingBlock.init_from_molecule(
                molecule=molecule,
                functional_groups=[
                    NPyridineFactory(), NTriazoleFactory()
                ]
            )

            # Update stk_mol to conformer geometry.
            new_mol = update_from_rdkit_conf(
                stk_mol=new_mol,
                rdk_mol=confs,
                conf_id=cid
            )

            angle = calculate_N_COM_N_angle(new_mol)

            if angle < min_angle:
                min_cid = cid
                min_angle = angle
                molecule = update_from_rdkit_conf(
                    stk_mol=molecule,
                    rdk_mol=confs,
                    conf_id=min_cid
                )

        xtb_opt = stko.XTB(
            xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
            output_dir=f'{name}',
            gfn_version=2,
            num_cores=6,
            opt_level='verytight',
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True
        )
        molecule = xtb_opt.optimize(mol=molecule)

    return molecule
