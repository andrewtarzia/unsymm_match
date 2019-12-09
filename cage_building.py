#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for building cages.

Author: Andrew Tarzia

Date Created: 5 Nov 2019

"""

from os.path import exists
from rdkit.Chem import AllChem as rdkit

import stk
import atools


def build_metal():
    m = rdkit.MolFromSmiles('[Pd+2]')
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
    metal_coord_info = {
        0: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        1: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        2: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
        3: {
            'atom_ids': [0],
            'bonder_ids': [0],
            'deleter_ids': [None]
        },
    }
    metal = stk.assign_metal_fgs(
        building_block=metal,
        coordination_info=metal_coord_info
    )
    return metal


def build_N_atom():
    m = rdkit.MolFromSmiles('N')
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    n_atom = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=['metal_bound_N'],
    )
    return n_atom


def build_metal_centre():
    """
    Build Pd2+ square planar metal centre.

    Parameters
    ----------

    Returns
    -------

    """

    metal = build_metal()
    n_atom = build_N_atom()
    sqpl = stk.metal_centre.SquarePlanar()
    complex = stk.ConstructedMolecule(
        building_blocks=[metal, n_atom],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            n_atom: sqpl.vertices[1:]
        }
    )
    complex = stk.BuildingBlock.init_from_molecule(
        complex,
        functional_groups=['metal_bound_N']
    )
    return complex


def optimize_cage(cage, cage_name):

    cage = atools.MOC_rdkit_opt(cage, cage_name, do_long=False)
    cage.write(f'{cage_name}_rdk.mol')
    cage.write(f'{cage_name}_rdk.xyz')

    cage = atools.MOC_uff_opt(cage, cage_name, metal_FFs={46: 'Pd4+2'})
    cage.write(f'{cage_name}_uff4mof.mol')
    cage.write(f'{cage_name}_uff4mof.xyz')

    cage = atools.MOC_MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=2,
        timestep='0.25',
        equib='0.5',
        production='0.5',
        metal_FFs={46: 'Pd4+2'},
        opt_conf=False,
        save_conf=False
    )

    cage = atools.MOC_MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=100,
        timestep='0.75',
        equib='0.5',
        production='100.0',
        metal_FFs={46: 'Pd4+2'},
        opt_conf=False,
        save_conf=True
    )

    cage.write(f'{cage_name}_MDconf.mol')
    cage.write(f'{cage_name}_MDconf.xyz')

    atools.MOC_xtb_conformers(
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
        solvent=('dmso', 'verytight')
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')

    atools.MOC_xtb_opt(
        cage,
        cage_name,
        nc=6,
        free_e=0,
        charge=4,
        opt_level='extreme',
        etemp=300,
        solvent=('dmso', 'verytight')
    )
    cage.write(f'{cage_name}_optc.mol')
    cage.write(f'{cage_name}_optc.xyz')
    cage.dump(f'{cage_name}_optc.json')


def build_cage_isomers(name, ligand, complex):
    """
    Build all four cage isomers.

    Parameters
    ----------

    Returns
    -------

    """

    cage_isomers = {}

    topologies = {
        'A': stk.cage.M2L4_Lantern(vertex_alignments={
            2: 0,
            3: 0,
            4: 0,
            5: 0
        }),
        'B': stk.cage.M2L4_Lantern(vertex_alignments={
            2: 1,
            3: 0,
            4: 0,
            5: 0
        }),
        'C': stk.cage.M2L4_Lantern(vertex_alignments={
            2: 1,
            3: 1,
            4: 0,
            5: 0
        }),
        'D': stk.cage.M2L4_Lantern(vertex_alignments={
            2: 1,
            3: 0,
            4: 1,
            5: 0
        }),
    }

    for top in topologies:
        topology = topologies[top]
        name_ = f'{name}_{top}'
        json_file = f'{name_}_optc.json'
        if exists(json_file):
            cage = stk.ConstructedMolecule.load(
                json_file
            )
        else:
            print(f'building {name_}')
            cage = stk.ConstructedMolecule(
                building_blocks=[complex, ligand],
                topology_graph=topology,
                building_block_vertices={
                    complex: topology.vertices[:2],
                    ligand: topology.vertices[2:],
                }
            )
            cage.write(f'{name_}_unopt.mol')
            cage.write(f'{name_}_unopt.xyz')
            cage.dump(f'{name_}_unopt.json')
            optimize_cage(cage, name_)

        cage_isomers[top] = cage

    return cage_isomers
