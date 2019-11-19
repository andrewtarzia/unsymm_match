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

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='stochastic',
        temperature='700',
        N=50,
        timestep='0.25',
        equib='0.5',
        production='5.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_MD.mol')
    cage.write(f'{cage_name}_MD.xyz')
    cage.dump(f'{cage_name}_MD.json')

    cage = uff_opt(cage, cage_name)

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='tight')


def optimize_cage_rdkit(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='tight')


def optimize_cage_rdkitMD(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='stochastic',
        temperature='700',
        N=50,
        timestep='0.25',
        equib='0.5',
        production='5.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='tight')


def opt_test1(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=True)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='700',
        N=100,
        timestep='0.25',
        equib='0.5',
        production='20.0',
        opt_conf=False
    )

    cage = uff_opt(cage, cage_name)

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='vtight')


def opt_test2(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='stochastic',
        temperature='700',
        N=100,
        timestep='0.25',
        equib='0.5',
        production='20.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='tight')


def opt_test3(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=10,
        timestep='0.25',
        equib='0.5',
        production='20.0',
        opt_conf=True
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='tight')


def opt_test4(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=200,
        timestep='0.25',
        equib='0.5',
        production='20.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='vtight')


def opt_test5(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=True)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=500,
        timestep='0.25',
        equib='0.5',
        production='30.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='extreme')


def opt_test6(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=True)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=10,
        timestep='0.25',
        equib='0.5',
        production='5.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_MD.mol')
    cage.write(f'{cage_name}_MD.xyz')
    cage.dump(f'{cage_name}_MD.json')

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=500,
        timestep='0.75',
        equib='0.5',
        production='100.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='extreme')


def opt_test7(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=10,
        timestep='0.25',
        equib='0.5',
        production='1.0',
        opt_conf=False
    )

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=1000,
        timestep='0.75',
        equib='0.5',
        production='100.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='extreme')


def opt_test8(cage, cage_name):

    cage = rdkit_opt(cage, cage_name, do_long=False)

    cage = uff_opt(cage, cage_name)

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='1000',
        N=10,
        timestep='0.25',
        equib='0.5',
        production='1.0',
        opt_conf=False
    )

    cage = MD_opt(
        cage,
        cage_name,
        integrator='leapfrog verlet',
        temperature='500',
        N=1000,
        timestep='0.75',
        equib='0.5',
        production='100.0',
        opt_conf=False
    )

    cage.write(f'{cage_name}_prextb.mol')
    cage.write(f'{cage_name}_prextb.xyz')
    cage.dump(f'{cage_name}_prextb.json')

    xtb_opt(cage, cage_name, opt_level='extreme')


def rdkit_opt(cage, cage_name, do_long):
    print('doing optimisation')
    optimizer = stk.MetalOptimizer(
        metal_binder_distance=1.9,
        metal_binder_fc=1.0e2,
        binder_ligand_fc=0.0,
        ignore_vdw=False,
        rel_distance=None,
        res_steps=50,
        restrict_bonds=True,
        restrict_angles=True,
        restrict_orientation=True,
        max_iterations=40,
        do_long_opt=do_long
    )

    optimizer.optimize(mol=cage)
    cage.write(f'{cage_name}_rdk.mol')
    cage.write(f'{cage_name}_rdk.xyz')
    cage.dump(f'{cage_name}_rdk.json')

    return cage


def uff_opt(cage, cage_name):

    metal_FFs = {46: 'Pd4+2'}

    print('doing UFF4MOF optimisation')
    gulp_opt = stk.GulpMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_FFs,
        output_dir=f'cage_opt_{cage_name}_uff1'
    )
    gulp_opt.assign_FF(cage)
    gulp_opt.optimize(mol=cage)
    cage.write(f'{cage_name}_uff4mof.mol')
    cage.write(f'{cage_name}_uff4mof.xyz')
    cage.dump(f'{cage_name}_uff4mof.json')

    return cage


def MD_opt(
    cage,
    cage_name,
    integrator,
    temperature,
    N,
    timestep,
    equib,
    production,
    opt_conf
):
    metal_FFs = {46: 'Pd4+2'}
    print('doing UFF4MOF MD')
    gulp_MD = stk.GulpMDMetalOptimizer(
        gulp_path='/home/atarzia/software/gulp-5.1/Src/gulp/gulp',
        metal_FF=metal_FFs,
        output_dir=f'cage_opt_{cage_name}_MD2',
        integrator=integrator,
        ensemble='nvt',
        temperature=temperature,
        equilbration=equib,
        production=production,
        timestep=timestep,
        N_conformers=N,
        opt_conformers=opt_conf
    )
    gulp_MD.assign_FF(cage)
    gulp_MD.optimize(cage)

    return cage


def xtb_opt(cage, cage_name, opt_level):
    print('doing XTB optimisation')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'cage_opt_{cage_name}_xtb',
        gfn_version=2,
        num_cores=6,
        opt_level=opt_level,
        charge=4,
        num_unpaired_electrons=0,
        max_runs=1,
        electronic_temperature=1000,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(mol=cage)
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
            # optimize_cage(cage, name_)
            # optimize_cage_rdkit(cage, name_)
            # optimize_cage_rdkitMD(cage, name_)
            # opt_test3(cage, name_)
            opt_test8(cage, name_)

        cage_isomers[top] = cage

    return cage_isomers
