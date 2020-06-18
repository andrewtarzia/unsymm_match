#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for building cages.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

from os.path import exists
import numpy as np

import stk
import atools


def optimize_cage(cage, cage_name):

    collapser_file = f'{cage_name}_coll.mol'
    uff4mof_file = f'{cage_name}_uff4mof.mol'
    mdconf_file = f'{cage_name}_MDconf.mol'
    prextb_file = f'{cage_name}_prextb.mol'
    opt_file = f'{cage_name}_optc.mol'

    if exists(collapser_file):
        cage = cage.with_structure_from_file(collapser_file)
    else:
        step_size = 0.05
        distance_cut = 2.0
        scale_steps = True
        cage = atools.MOC_collapse(
            cage,
            cage_name,
            step_size=step_size,
            distance_cut=distance_cut,
            scale_steps=scale_steps
        )
        cage.write(collapser_file)

    if exists(uff4mof_file):
        cage = cage.with_structure_from_file(uff4mof_file)
    else:
        cage = atools.MOC_uff_opt(
            cage,
            cage_name,
            metal_FFs={46: 'Pd4+2'}
        )
        cage.write(uff4mof_file)

    if exists(mdconf_file):
        cage = cage.with_structure_from_file(mdconf_file)
    else:
        cage = atools.MOC_MD_opt(
            cage,
            cage_name,
            integrator='leapfrog verlet',
            temperature=1000,
            N=2,
            timestep=0.25,
            equib=0.5,
            production=0.5,
            metal_FFs={46: 'Pd4+2'},
            opt_conf=False,
            save_conf=False
        )

        cage = atools.MOC_MD_opt(
            cage,
            cage_name,
            integrator='leapfrog verlet',
            temperature=1000,
            N=100,
            timestep=0.75,
            equib=0.5,
            production=100.0,
            metal_FFs={46: 'Pd4+2'},
            opt_conf=False,
            save_conf=True
        )

        cage.write(mdconf_file)

    if exists(prextb_file):
        cage = cage.with_structure_from_file(prextb_file)
    else:
        cage = atools.MOC_xtb_conformers(
            cage,
            cage_name,
            opt=True,
            opt_level='normal',
            nc=6,
            free_e=0,
            charge=4,
            etemp=300,
            conformer_dir=f'cage_opt_{cage_name}_MD',
            output_dir=f'cage_opt_{cage_name}_xtb_conf',
            solvent=('dmso', 'verytight'),
            handle_failure=True
        )
        cage.write(prextb_file)

    cage = atools.MOC_xtb_opt(
        cage,
        cage_name,
        nc=6,
        free_e=0,
        charge=4,
        opt_level='extreme',
        etemp=300,
        solvent=('dmso', 'verytight')
    )
    cage.write(opt_file)

    return cage


def build_cage_isomers(name, ligand):
    """
    Build all four cage isomers.

    Parameters
    ----------

    Returns
    -------

    """

    cage_isomers = {}

    v_alignments = {
        'A': {2: 0, 3: 0, 4: 0, 5: 0},
        'B': {2: 1, 3: 0, 4: 0, 5: 0},
        'C': {2: 1, 3: 1, 4: 0, 5: 0},
        'D': {2: 1, 3: 0, 4: 1, 5: 0},
    }

    complex = stk.BuildingBlock(
        smiles='[Pd+2]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2))
            for i in range(4)
        ),
        position_matrix=np.array([[0, 0, 0]]),
    )

    for v_a in v_alignments:
        v_align = v_alignments[v_a]
        name_ = f'{name}_{v_a}'
        opt_file = f'{name_}_optc.mol'

        # Build cage.
        cage = stk.ConstructedMolecule(
            stk.cage.M2L4Lantern(
                building_blocks={
                    complex: (0, 1),
                    ligand: (2, 3, 4, 5)
                },
                vertex_alignments=v_align,
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                atools.AromaticCNC,
                                stk.SingleAtom
                            }): 9,
                            frozenset({
                                atools.AromaticCNN,
                                stk.SingleAtom
                            }): 9,
                        }
                    )
                )
            )
        )

        if exists(opt_file):
            print(
                'Using non constructed molecule here. fix this in '
                'future'
            )
            # cage = cage.with_structure_from_file(opt_file)
            cage = stk.BuildingBlock.init_from_file(opt_file)
        else:
            print(f'optimizing {name_}')
            cage.write(f'{name_}_unopt.mol')
            cage = optimize_cage(cage, name_)

        cage_isomers[v_a] = cage

    return cage_isomers
