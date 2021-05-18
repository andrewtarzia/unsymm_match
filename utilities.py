#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utilities for screening.

Author: Andrew Tarzia

Date Created: 18 May 2021

"""

import numpy as np
import networkx as nx
from itertools import combinations
import glob
from scipy.spatial.distance import euclidean
from mendeleev import element
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem as rdkit
from rdkit.Chem import Draw
import os
import json
import pymatgen as pmg
from pymatgen.analysis.local_env import (
    LocalStructOrderParams,
)
from pymatgen.core import Specie

import stk
import stko

from plotting import scatter_plot


def update_from_rdkit_conf(stk_mol, rdk_mol, conf_id):
    """
    Update the structure to match `conf_id` of `mol`.

    Parameters
    ----------
    struct : :class:`stk.Molecule`
        The molecule whoce coordinates are to be updated.

    mol : :class:`rdkit.Mol`
        The :mod:`rdkit` molecule to use for the structure update.

    conf_id : :class:`int`
        The conformer ID of the `mol` to update from.

    Returns
    -------
    :class:`.Molecule`
        The molecule.

    """

    pos_mat = rdk_mol.GetConformer(id=conf_id).GetPositions()
    return stk_mol.with_position_matrix(pos_mat)


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249
    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, normal=None):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    If normal is given, the angle polarity is determined using the
    cross product of the two vectors.

    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if normal is not None:
        # Get normal vector and cross product to determine sign.
        cross = np.cross(v1_u, v2_u)
        if np.dot(normal, cross) < 0:
            angle = -angle
    return angle


def get_atom_distance(molecule, atom1_id, atom2_id):
    """
    Return the distance between atom1 and atom2.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`

    atom1_id : :class:`int`
        The id of atom1.

    atom2_id : :class:`int`
        The id of atom2.

    Returns
    -------
    :class:`float`
        The euclidean distance between two atoms.

    """

    position_matrix = molecule.get_position_matrix()

    distance = euclidean(
        u=position_matrix[atom1_id],
        v=position_matrix[atom2_id]
    )

    return float(distance)


def calculate_NN_distance(bb):
    """
    Calculate the N-N distance of ditopic building block.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC and metals + AromaticCNN.

    Parameters
    ----------
    bb : :class:`stk.BuildingBlock`
        stk molecule to analyse.

    Returns
    -------
    NN_distance : :class:`float`
        Distance(s) between [angstrom] N atoms in functional groups.

    """

    fg_counts = 0
    N_positions = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC) or isinstance(fg, AromaticCNN):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            N_position, = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            N_positions.append(N_position)

    if fg_counts != 2:
        raise ValueError(
            f'{bb} does not have 2 AromaticCNC or AromaticCNN '
            'functional groups.'
        )

    # Calculate the angle between the two vectors.
    NN_distance = np.linalg.norm(N_positions[0] - N_positions[1])
    return NN_distance


def calculate_bite_angle(bb):
    """
    Calculate the bite angle of a ditopic building block.

    Here the bite angle is defined `visually` as in:
        https://doi.org/10.1016/j.ccr.2018.06.010

    It is calculated using the angles between binding vectors and the
    N to N vector.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC and metals + AromaticCNN.

    Parameters
    ----------
    bb : :class:`stk.BuildingBlock`
        Stk molecule to analyse.

    Returns
    -------
    bite_angle : :class:`float`
        Angle between two bonding vectors of molecule.

    """

    fg_counts = 0
    fg_vectors = []
    N_positions = []

    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC) or isinstance(fg, AromaticCNN):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position.
            N_position, = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            N_positions.append(N_position)
            # Get centroid of neighbouring C atom positions.
            if isinstance(fg, AromaticCNC):
                CC_MP = bb.get_centroid(
                    atom_ids=(
                        fg.get_carbon1().get_id(),
                        fg.get_carbon2().get_id()
                    )
                )
            elif isinstance(fg, AromaticCNN):
                CC_MP = bb.get_centroid(
                    atom_ids=(
                        fg.get_carbon().get_id(),
                        fg.get_nitrogen2().get_id()
                    )
                )
            # Get vector between COM and N position.
            v = N_position - CC_MP
            fg_vectors.append(v)

    if fg_counts != 2:
        raise ValueError(
            f'{bb} does not have exactly 2 AromaticCNC or AromaticCNN '
            'functional groups.'
        )

    # Get N to N vector.
    NN_vec = N_positions[1] - N_positions[0]

    # Calculate the angle between the two vectors.
    angle_1 = np.degrees(angle_between(
        fg_vectors[0], NN_vec
    ))
    angle_2 = np.degrees(angle_between(
        fg_vectors[1], -NN_vec
    ))
    bite_angle = (angle_1 - 90) + (angle_2 - 90)
    return bite_angle


def get_center_of_mass(molecule, atom_ids=None):
    """
    Return the centre of mass.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`

    atom_ids : :class:`iterable` of :class:`int`, optional
        The ids of atoms which should be used to calculate the
        center of mass. If ``None``, then all atoms will be used.

    Returns
    -------
    :class:`numpy.ndarray`
        The coordinates of the center of mass.

    References
    ----------
    https://en.wikipedia.org/wiki/Center_of_mass

    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(atom_ids, (list, tuple)):
        # Iterable gets used twice, once in get_atom_positions
        # and once in zip.
        atom_ids = list(atom_ids)

    center = 0
    total_mass = 0.
    coords = molecule.get_atomic_positions(atom_ids)
    atoms = molecule.get_atoms(atom_ids)
    for atom, coord in zip(atoms, coords):
        mass = element(atom.__class__.__name__).atomic_weight
        total_mass += mass
        center += mass*coord
    return np.divide(center, total_mass)


def calculate_N_COM_N_angle(bb):
    """
    Calculate the N-COM-N angle of a ditopic building block.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC and metals + AromaticCNN.

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
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC) or isinstance(fg, AromaticCNN):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            N_position, = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            fg_positions.append(N_position)

    if fg_counts != 2:
        raise ValueError(
            f'{bb} does not have 2 AromaticCNC or AromaticCNN '
            'functional groups.'
        )

    # Get building block COM.
    COM_position = get_center_of_mass(bb)

    # Get vectors.
    fg_vectors = [i-COM_position for i in fg_positions]

    # Calculate the angle between the two vectors.
    angle = np.degrees(angle_between(*fg_vectors))
    return angle


def draw_and_save_grid(
    mol_list,
    names,
    subImgSize,
    mol_per_row,
    filename
):
    """
    Draw RDKit molecules and save SVG.

    """
    img = Draw.MolsToGridImage(
        mol_list,
        molsPerRow=mol_per_row,
        subImgSize=subImgSize,
        legends=names,
        useSVG=True
    )
    save_svg(
        filename=filename,
        string=img
    )


def draw_molecules(
    ligands,
    energy_preferences,
    plane_devs,
    sqpl_ops,
):
    """
    Draw molecules as grids with scores.

    """

    mol_list = []
    name_list = []
    count = 0
    for i, lig in enumerate(ligands):
        MOL = rdkit.MolFromMolFile(
            f'{lig}_opt.mol'
        )
        MOL.RemoveAllConformers()
        mol_list.append(MOL)
        name_list.append(
            f'{lig}: OP={round(sqpl_ops[i], 2)}'
            f',PD={round(plane_devs[i], 2)}\n'
            f'E={round(energy_preferences[i], 2)}'
        )
        count += 1

    # Sort by energy preferences.
    mol_list = [
        x for _, x in sorted(zip(energy_preferences, mol_list))
    ]
    name_list = [
        x for _, x in sorted(zip(energy_preferences, name_list))
    ]

    if not os.path.exists('molecules_scoring'):
        os.mkdir('molecules_scoring')

    # Save figure of desired molecules.
    mol_list2grid(
        molecules=mol_list,
        names=name_list,
        filename='molecules_scoring/molecules_scoring',
        mol_per_row=4,
        maxrows=3,
        subImgSize=(300, 300)
    )


def save_svg(filename, string):
    """
    Save svg text to a file.

    """

    with open(filename, 'w') as f:
        f.write(string)


def mol_list2grid(
    molecules,
    filename,
    mol_per_row,
    maxrows,
    subImgSize=(200, 200),
    names=None
):
    """
    Produce a grid of molecules in mol_list.

    molecules (list) - list of molecule SMILEs

    """

    if len(molecules) > mol_per_row * maxrows:
        # have to make multiple images
        new_mol_list = []
        new_names = []
        count = 1
        for i, mol in enumerate(molecules):
            new_mol_list.append(mol)
            if names is None:
                new_names = None
            else:
                new_names.append(names[i])
            # make image
            chk1 = len(new_mol_list) == mol_per_row * maxrows
            chk2 = i == len(molecules)-1
            if chk1 or chk2:
                draw_and_save_grid(
                    mol_list=new_mol_list,
                    mol_per_row=mol_per_row,
                    subImgSize=subImgSize,
                    names=new_names,
                    filename=f'{filename}_{count}.svg'
                )
                # img.save(filename + '_' + str(count) + '.png')
                new_mol_list = []
                new_names = []
                count += 1
    else:
        draw_and_save_grid(
            mol_list=molecules,
            mol_per_row=mol_per_row,
            subImgSize=subImgSize,
            names=names,
            filename=f'{filename}.svg'
        )


def MOC_collapse(
    cage,
    cage_name,
    step_size,
    distance_cut,
    scale_steps
):
    """
    Perform Collapser optimisation of MOC.

    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.

    cage_name : :class:`str`
        Name of cage.

    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.

    """

    print(f'..........doing collapser optimisation of {cage_name}')
    output_dir = f'cage_opt_{cage_name}_coll'
    optimizer = stko.Collapser(
        output_dir=output_dir,
        step_size=step_size,
        distance_cut=distance_cut,
        scale_steps=scale_steps,
    )
    cage = optimizer.optimize(mol=cage)

    return cage


def MOC_uff_opt(
    cage,
    cage_name,
    metal_FFs,
    metal_ligand_bond_order='half',
    CG=False,
    maxcyc=1000,
    gulp_exec=None,
):
    """
    Perform UFF4MOF optimisation of MOC.

    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.

    cage_name : :class:`str`
        Name of cage.

    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.

    """

    if gulp_exec is None:
        gulp_exec = '/home/atarzia/software/gulp-5.1/Src/gulp/gulp'

    output_dir = (
        f'cage_opt_{cage_name}_uff' if CG is False
        else f'cage_opt_{cage_name}_uffCG'
    )

    print(f'..........doing UFF4MOF optimisation of {cage_name}')
    print(f'Conjugate Gradient: {CG}, Max steps: {maxcyc}')
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path=gulp_exec,
        maxcyc=maxcyc,
        metal_FF=metal_FFs,
        metal_ligand_bond_order=metal_ligand_bond_order,
        output_dir=output_dir,
        conjugate_gradient=CG
    )
    gulp_opt.assign_FF(cage)
    cage = gulp_opt.optimize(mol=cage)

    return cage


def MOC_MD_opt(
    cage,
    cage_name,
    integrator,
    temperature,
    N,
    timestep,
    equib,
    production,
    opt_conf,
    metal_FFs,
    metal_ligand_bond_order='half',
    save_conf=False,
    gulp_exec=None,
):
    """
    Perform UFF4MOF molecular dynamics of MOC.

    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.

    cage_name : :class:`str`
        Name of cage.

    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.

    """

    if gulp_exec is None:
        gulp_exec = '/home/atarzia/software/gulp-5.1/Src/gulp/gulp'

    print(f'..........doing UFF4MOF MD of {cage_name}')
    gulp_MD = stko.GulpUFFMDOptimizer(
        gulp_path=gulp_exec,
        metal_FF=metal_FFs,
        metal_ligand_bond_order=metal_ligand_bond_order,
        output_dir=f'cage_opt_{cage_name}_MD',
        integrator=integrator,
        ensemble='nvt',
        temperature=temperature,
        equilbration=equib,
        production=production,
        timestep=timestep,
        N_conformers=N,
        opt_conformers=opt_conf,
        save_conformers=save_conf
    )
    gulp_MD.assign_FF(cage)
    cage = gulp_MD.optimize(cage)

    return cage


def MOC_xtb_conformers(
    cage,
    cage_name,
    etemp,
    output_dir,
    conformer_dir,
    nc,
    free_e,
    charge,
    gfn_exec=None,
    opt=False,
    opt_level=None,
    solvent=None,
    handle_failure=False
):
    """
    Perform GFN2-xTB conformer scan of MOC.

    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.

    cage_name : :class:`str`
        Name of cage.

    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.

    """

    if gfn_exec is None:
        gfn_exec = '/home/atarzia/software/xtb-190806/bin/xtb'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print(
        f'..........doing XTB conformer sorting by '
        f'energy of {cage_name}'
    )
    conformers = glob.glob(f'{conformer_dir}/conf_*.xyz')
    ids = []
    energies = []
    min_energy = 10E20
    for file in sorted(conformers):
        id = file.replace('.xyz', '').split('_')[-1]
        cage = cage.with_structure_from_file(file)
        opt_failed = False
        if opt:
            print(f'optimising conformer {id}')
            xtb_opt = stko.XTB(
                xtb_path=gfn_exec,
                output_dir=f'opt_{cage_name}_{id}',
                gfn_version=2,
                num_cores=nc,
                opt_level=opt_level,
                charge=charge,
                num_unpaired_electrons=free_e,
                max_runs=1,
                electronic_temperature=etemp,
                calculate_hessian=False,
                unlimited_memory=True,
                solvent=solvent_str,
                solvent_grid=solvent_grid
            )
            try:
                cage = xtb_opt.optimize(mol=cage)
                cage.write(os.path.join(
                    f'{output_dir}',
                    f'conf_{id}_opt.xyz',
                ))
            except stko.XTBConvergenceError:
                if handle_failure:
                    opt_failed = True
                    print(f'optimising conformer {id}: FAILED')
                else:
                    raise stko.XTBConvergenceError()

        print(f'..........calculating energy of {id} of {cage_name}')
        # Extract energy.
        xtb_energy = stko.XTBEnergy(
            xtb_path=gfn_exec,
            output_dir=f'ey_{cage_name}_{id}',
            num_cores=nc,
            charge=charge,
            num_unpaired_electrons=free_e,
            electronic_temperature=etemp,
            unlimited_memory=True,
            solvent=solvent_str,
            solvent_grid=solvent_grid
        )
        if handle_failure and opt_failed:
            energy = 10E24
        else:
            energy = xtb_energy.get_energy(cage)
        if energy < min_energy:
            min_energy_conformer = file
            min_energy = energy
        ids.append(id)
        energies.append(energy)

    print('done', min_energy, min_energy_conformer)
    cage = cage.with_structure_from_file(min_energy_conformer)

    energies = [(i-min(energies))*2625.5 for i in energies]
    fig, ax = scatter_plot(
        X=ids, Y=energies,
        xtitle='conformer id',
        ytitle='rel. energy [kJmol$^{-1}$]',
        xlim=(0, 201),
        ylim=(-5, 1000)
    )

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, f'{cage_name}_conf_energies.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    return cage


def MOC_xtb_opt(
    cage,
    cage_name,
    nc,
    opt_level,
    etemp,
    charge,
    free_e,
    gfn_exec=None,
    solvent=None
):
    """
    Perform GFN2-xTB optimisation of MOC.

    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.

    cage_name : :class:`str`
        Name of cage.

    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.

    """

    if gfn_exec is None:
        gfn_exec = '/home/atarzia/software/xtb-190806/bin/xtb'

    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print(f'..........doing XTB optimisation of {cage_name}')
    xtb_opt = stko.XTB(
        xtb_path=gfn_exec,
        output_dir=f'cage_opt_{cage_name}_xtb',
        gfn_version=2,
        num_cores=nc,
        opt_level=opt_level,
        charge=charge,
        num_unpaired_electrons=free_e,
        max_runs=1,
        electronic_temperature=etemp,
        calculate_hessian=False,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    cage = xtb_opt.optimize(mol=cage)

    return cage


class AromaticCNCFactory(stk.FunctionalGroupFactory):
    """
    A subclass of stk.SmartsFunctionalGroupFactory.

    """

    def __init__(self, bonders=(1, ), deleters=()):
        """
        Initialise :class:`.AromaticCNCFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        generic_functional_groups = stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#6]',
            bonders=self._bonders,
            deleters=self._deleters
        ).get_functional_groups(molecule)
        for fg in generic_functional_groups:
            atom_ids = (i.get_id() for i in fg.get_atoms())
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield AromaticCNC(
                carbon1=atoms[0],
                nitrogen=atoms[1],
                carbon2=atoms[2],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )


class AromaticCNC(stk.GenericFunctionalGroup):
    """
    Represents an N atom in pyridine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon][nitrogen][carbon]``.

    """

    def __init__(self, carbon1, nitrogen, carbon2, bonders, deleters):
        """
        Initialize a :class:`.Alcohol` instance.

        Parameters
        ----------
        carbon1 : :class:`.C`
            The first carbon atom.

        nitrogen : :class:`.N`
            The nitrogen atom.

        carbon2 : :class:`.C`
            The second carbon atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon1 = carbon1
        self._nitrogen = nitrogen
        self._carbon2 = carbon2
        atoms = (carbon1, nitrogen, carbon2)
        super().__init__(atoms, bonders, deleters)

    def get_carbon1(self):
        """
        Get the first carbon atom.

        Returns
        -------
        :class:`.C`
            The first carbon atom.

        """

        return self._carbon1

    def get_carbon2(self):
        """
        Get the second carbon atom.

        Returns
        -------
        :class:`.C`
            The second carbon atom.

        """

        return self._carbon2

    def get_nitrogen(self):
        """
        Get the nitrogen atom.

        Returns
        -------
        :class:`.N`
            The nitrogen atom.

        """

        return self._nitrogen

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._nitrogen = self._nitrogen
        clone._carbon2 = self._carbon2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._nitrogen}, {self._carbon2}, '
            f'bonders={self._bonders})'
        )


class AromaticCNNFactory(stk.FunctionalGroupFactory):
    """
    A subclass of stk.SmartsFunctionalGroupFactory.

    """

    def __init__(self, bonders=(1, ), deleters=()):
        """
        Initialise :class:`.AromaticCNNFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        generic_functional_groups = stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#7X2]',
            bonders=self._bonders,
            deleters=self._deleters
        ).get_functional_groups(molecule)
        for fg in generic_functional_groups:
            atom_ids = (i.get_id() for i in fg.get_atoms())
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield AromaticCNN(
                carbon=atoms[0],
                nitrogen=atoms[1],
                nitrogen2=atoms[2],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )


class AromaticCNN(stk.GenericFunctionalGroup):
    """
    Represents an N atom in pyridine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon][nitrogen][nitrogen]``.

    """

    def __init__(
        self,
        carbon,
        nitrogen,
        nitrogen2,
        bonders,
        deleters
    ):
        """
        Initialize a :class:`.Alcohol` instance.

        Parameters
        ----------
        carbon : :class:`.C`
            The carbon atom.

        nitrogen : :class:`.N`
            The first and bonding (default) nitrogen atom.

        nitrogen2 : :class:`.C`
            The second nitrogen atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon = carbon
        self._nitrogen = nitrogen
        self._nitrogen2 = nitrogen2
        atoms = (carbon, nitrogen, nitrogen2)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        """
        Get the carbon atom.

        Returns
        -------
        :class:`.C`
            The carbon atom.

        """

        return self._carbon

    def get_nitrogen2(self):
        """
        Get the second nitrogen atom.

        Returns
        -------
        :class:`.N`
            The second nitrogen atom.

        """

        return self._nitrogen2

    def get_nitrogen(self):
        """
        Get the first nitrogen atom.

        Returns
        -------
        :class:`.N`
            The first nitrogen atom.

        """

        return self._nitrogen

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._nitrogen = self._nitrogen
        clone._nitrogen2 = self._nitrogen2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._nitrogen2 = atom_map.get(
            self._nitrogen2.get_id(),
            self._nitrogen2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._nitrogen}, {self._nitrogen2}, '
            f'bonders={self._bonders})'
        )


def get_dihedral(pt1, pt2, pt3, pt4):
    """
    Calculate the dihedral between four points.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product

    Output in range (-pi to pi).

    From: https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-
    coordinates-in-python
    (new_dihedral(p))

    """

    p0 = np.asarray(pt1)
    p1 = np.asarray(pt2)
    p2 = np.asarray(pt3)
    p3 = np.asarray(pt4)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def get_stk_bond_angle(mol, atom1_id, atom2_id, atom3_id):
    # TODO: Fill in the doc string for this including defintiions.
    atom1_pos, = mol.get_atomic_positions(atom_ids=atom1_id)
    atom2_pos, = mol.get_atomic_positions(atom_ids=atom2_id)
    atom3_pos, = mol.get_atomic_positions(atom_ids=atom3_id)
    v1 = atom1_pos - atom2_pos
    v2 = atom3_pos - atom2_pos
    return stk.vector_angle(v1, v2)


def shortest_distance_to_plane(plane, point):
    """
    Calculate the perpendicular distance beween a point and a plane.

    """

    top = abs(
        plane[0]*point[0] + plane[1]*point[1] +
        plane[2]*point[2] - plane[3]
    )
    bottom = np.sqrt(plane[0]**2 + plane[1]**2 + plane[2]**2)
    distance = top / bottom
    return distance


def get_square_planar_distortion(mol, metal, bonder):
    """
    Calculate measures of distortion of a square planer metal.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal : :class:`int`
        Element number of metal atom.

    bonder : :class:`int`
        Element number of atoms bonded to metal.

    Returns
    -------
    results : :class:`dict`
        Dictionary containing 'bond_lengths', 'angles', 'torsions' and
        'plane_dev'.

    """

    results = {
        'bond_lengths': [],
        'angles': [],
        'torsions': [],
        'plane_dev': [],
        'plane_angle_avg': [],
        'plane_angle_std': []
    }

    # Find metal atoms.
    metal_atoms = []
    for atom in mol.get_atoms():
        if atom.get_atomic_number() == metal:
            metal_atoms.append(atom)

    # Find bonders.
    metal_bonds = []
    ids_to_metals = []
    for bond in mol.get_bonds():
        if bond.get_atom1() in metal_atoms:
            metal_bonds.append(bond)
            ids_to_metals.append(bond.get_atom2().get_id())
        elif bond.get_atom2() in metal_atoms:
            metal_bonds.append(bond)
            ids_to_metals.append(bond.get_atom1().get_id())

    # Calculate bond lengths.
    for bond in metal_bonds:
        results['bond_lengths'].append(
            get_atom_distance(
                molecule=mol,
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
            )
        )

    # Calculate bond angles.
    for bonds in combinations(metal_bonds, r=2):
        bond1, bond2 = bonds
        bond1_atoms = [bond1.get_atom1(), bond1.get_atom2()]
        bond2_atoms = [bond2.get_atom1(), bond2.get_atom2()]
        pres_atoms = list(set(bond1_atoms + bond2_atoms))
        # If there are more than 3 atoms, implies two
        # independant bonds.
        if len(pres_atoms) > 3:
            continue
        for atom in pres_atoms:
            if atom in bond1_atoms and atom in bond2_atoms:
                idx2 = atom.get_id()
            elif atom in bond1_atoms:
                idx1 = atom.get_id()
            elif atom in bond2_atoms:
                idx3 = atom.get_id()

        angle = np.degrees(get_stk_bond_angle(
            mol=mol,
            atom1_id=idx1,
            atom2_id=idx2,
            atom3_id=idx3,
        ))
        if angle < 120:
            results['angles'].append(angle)

    # Calculate torsion.
    for metal_atom in metal_atoms:
        torsion_ids = []
        for bond in metal_bonds:
            if metal_atom.get_id() == bond.get_atom1().get_id():
                torsion_ids.append(bond.get_atom2().get_id())
            elif metal_atom.get_id() == bond.get_atom2().get_id():
                torsion_ids.append(bond.get_atom1().get_id())
        atom_positions = [
            i for i in mol.get_atomic_positions(atom_ids=torsion_ids)
        ]
        results['torsions'].append(abs(get_dihedral(*atom_positions)))

    # Calculate deviation of metal from bonder plane.
    for metal_atom in metal_atoms:
        binder_atom_ids = [metal_atom.get_id()]
        for bond in metal_bonds:
            if metal_atom.get_id() == bond.get_atom1().get_id():
                binder_atom_ids.append(bond.get_atom2().get_id())
            elif metal_atom.get_id() == bond.get_atom2().get_id():
                binder_atom_ids.append(bond.get_atom1().get_id())
        centroid = mol.get_centroid(atom_ids=binder_atom_ids)
        normal = mol.get_plane_normal(atom_ids=binder_atom_ids)
        # Plane of equation ax + by + cz = d.
        binder_atom_plane = np.append(normal, np.sum(normal*centroid))
        # Define the plane deviation as the sum of the distance of all
        # atoms from the plane defined by all atoms.
        plane_dev = sum([
            shortest_distance_to_plane(
                binder_atom_plane,
                tuple(mol.get_atomic_positions(atom_ids=i), )[0]
            )
            for i in binder_atom_ids
        ])
        results['plane_dev'].append(plane_dev)

    # Calculate N-Pd bond vector and CNC plane angle.
    for metal_atom in metal_atoms:
        plane_angles = []
        # Get metal position.
        metal_position = mol.get_centroid(
            atom_ids=[metal_atom.get_id()]
        )
        # Iterate over metal bonds.
        for bond in metal_bonds:
            if metal_atom.get_id() == bond.get_atom1().get_id():
                N_atom_id = bond.get_atom2().get_id()
            elif metal_atom.get_id() == bond.get_atom2().get_id():
                N_atom_id = bond.get_atom1().get_id()
            else:
                continue
            # Get MN vector.
            N_position = mol.get_centroid(atom_ids=[N_atom_id])
            MN_vector = N_position - metal_position
            # Get CNC atom ids.
            CNC_atom_ids = [N_atom_id]
            for bond in mol.get_bonds():
                if metal_atom.get_id() in [
                    bond.get_atom1().get_id(),
                    bond.get_atom2().get_id()
                ]:
                    continue
                if bond.get_atom1().get_id() == N_atom_id:
                    CNC_atom_ids.append(bond.get_atom2().get_id())
                elif bond.get_atom2().get_id() == N_atom_id:
                    CNC_atom_ids.append(bond.get_atom1().get_id())

            # Get CNC plane.
            centroid = mol.get_centroid(atom_ids=CNC_atom_ids)
            CNC_plane_normal = mol.get_plane_normal(
                atom_ids=CNC_atom_ids
            )
            # Calculate angle between CNC plane and MN vector.
            pa = np.degrees(angle_between(MN_vector, CNC_plane_normal))
            plane_angles.append(pa)

        # Define the plane angle of a metal centre as the sum of all
        # plane angles of 4 coordinated atoms.
        plane_angle_avg = np.average([i for i in plane_angles])
        plane_angle_std = np.std([i for i in plane_angles])
        results['plane_angle_avg'].append(plane_angle_avg)
        results['plane_angle_std'].append(plane_angle_std)

    return results


def convert_stk_to_pymatgen(stk_mol):
    """
    Convert stk.Molecule to pymatgen.Molecule.

    Parameters
    ----------
    stk_mol : :class:`stk.Molecule`
        Stk molecule to convert.

    Returns
    -------
    pmg_mol : :class:`pymatgen.Molecule`
        Corresponding pymatgen Molecule.

    """
    stk_mol.write('temp.xyz')
    pmg_mol = pmg.Molecule.from_file('temp.xyz')
    os.system('rm temp.xyz')

    return pmg_mol


def calculate_sites_order_values(
    molecule,
    site_idxs,
    target_species_type=None,
    neigh_idxs=None
):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    molecule : :class:`pmg.Molecule` or :class:`pmg.Structure`
        Pymatgen (pmg) molecule/structure to analyse.

    site_idxs : :class:`list` of :class:`int`
        Atom ids of sites to calculate OP of.

    target_species_type : :class:`str`
        Target neighbour element to use in OP calculation.
        Defaults to :class:`NoneType` if no target species is known.

    neigh_idxs : :class:`list` of :class:`list` of :class:`int`
        Neighbours of each atom in site_idx. Ordering is important.
        Defaults to :class:`NoneType` for when using
        :class:`pmg.Structure` - i.e. a structure with a lattice.

    Returns
    -------
    results : :class:`dict`
        Dictionary of format
        site_idx: dict of order parameters
        {
            `oct`: :class:`float`,
            `sq_plan`: :class:`float`,
            `q2`: :class:`float`,
            `q4`: :class:`float`,
            `q6`: :class:`float`
        }.

    """

    results = {}

    if target_species_type is None:
        targ_species = None
    else:
        targ_species = Specie(target_species_type)

    # Define local order parameters class based on desired types.
    types = [
        'oct',  # Octahedra OP.
        'sq_plan',  # Square planar envs.
        'q2',  # l=2 Steinhardt OP.
        'q4',  # l=4 Steinhardt OP.
        'q6',  # l=6 Steinhardt OP.
    ]
    loc_ops = LocalStructOrderParams(
        types=types,
    )
    if neigh_idxs is None:
        for site in site_idxs:
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                target_spec=[targ_species]
            )
            results[site] = {i: j for i, j in zip(types, site_results)}
    else:
        for site, neigh in zip(site_idxs, neigh_idxs):
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                indices_neighs=neigh,
                target_spec=targ_species
            )
            results[site] = {i: j for i, j in zip(types, site_results)}

    return results


def get_order_values(mol, metal, per_site=False):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal : :class:`int`
        Element number of metal atom.

    per_site : :class:`bool`
        Defaults to False. True if the OPs for each site are desired.

    Returns
    -------
    results : :class:`dict`
        Dictionary of order parameter max/mins/averages if `per_site`
        is False.

    """

    pmg_mol = convert_stk_to_pymatgen(stk_mol=mol)
    # Get sites of interest and their neighbours.
    sites = []
    neighs = []
    for atom in mol.get_atoms():
        if atom.get_atomic_number() == metal:
            sites.append(atom.get_id())
            bonds = [
                i
                for i in mol.get_bonds()
                if i.get_atom1().get_id() == atom.get_id()
                or i.get_atom2().get_id() == atom.get_id()
            ]
            a_neigh = []
            for b in bonds:
                if b.get_atom1().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom2().get_id())
                elif b.get_atom2().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom1().get_id())
            neighs.append(a_neigh)

    order_values = calculate_sites_order_values(
        molecule=pmg_mol,
        site_idxs=sites,
        neigh_idxs=neighs
    )

    if per_site:
        results = order_values
        return results
    else:
        # Get max, mins and averages of all OPs for the whole molecule.
        OPs = [order_values[i].keys() for i in order_values][0]
        OP_lists = {}
        for OP in OPs:
            OP_lists[OP] = [order_values[i][OP] for i in order_values]

        results = {
            # OP: (min, max, avg)
            i: {
                'min': min(OP_lists[i]),
                'max': max(OP_lists[i]),
                'avg': np.average(OP_lists[i])
            }
            for i in OP_lists
        }

        return results


def read_gfnx2xtb_eyfile(file):
    """
    Read the energy (kJ/mol from GFN2-xTB) from a .ey file.

    """

    with open(file, 'r') as f:
        lines = f.readlines()
        ey = float(lines[0].rstrip())

    return ey*2625.5


def calculate_energy(
    name,
    mol,
    ey_file,
    gfn_exec=None,
    charge=0,
    no_unpaired_e=0,
    solvent=None
):
    """
    Calculate GFN-xTB energy of molecule.

    """

    if gfn_exec is None:
        gfn_exec = '/home/atarzia/software/xtb-190806/bin/xtb'

    print(f'....getting energy of {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    xtb_energy = stko.XTBEnergy(
        xtb_path=gfn_exec,
        output_dir=f'{name}_ey',
        num_cores=6,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        electronic_temperature=300,
        unlimited_memory=True,
        calculate_free_energy=False,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    energy = xtb_energy.get_energy(mol)

    with open(ey_file, 'w') as f:
        f.write(str(energy))


def calculate_ligand_SE(
    org_ligs,
    smiles_keys,
    output_json,
    file_prefix=None
):
    """
    Calculate the strain energy of each ligand in the cage.

    Parameters
    ----------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    output_json : :class:`str`
        File name to save output to to avoid reruns.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns
    -------
    strain_energies : :class:`dict`
        Strain energies for each ligand.

    """

    # Check if output file exists.
    if not os.path.exists(output_json):
        strain_energies = {}
        # Iterate over ligands.
        for lig in org_ligs:
            stk_lig = org_ligs[lig]
            ey_file = lig.replace('mol', 'ey')
            smiles_key = stk.Smiles().get_key(stk_lig)
            idx = smiles_keys[smiles_key]
            sgt = str(stk_lig.get_num_atoms())
            # Get optimized ligand name that excludes any cage
            # information.
            if file_prefix is None:
                filename_ = f'organic_linker_s{sgt}_{idx}_opt.mol'
                opt_lig_ey = f'organic_linker_s{sgt}_{idx}_opt.ey'
                opt_lig_n = f'organic_linker_s{sgt}_{idx}_opt'
            else:
                filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'
                opt_lig_ey = f'{file_prefix}{sgt}_{idx}_opt.ey'
                opt_lig_n = f'{file_prefix}{sgt}_{idx}_opt'

            # Calculate energy of extracted ligand.
            if not os.path.exists(ey_file):
                calculate_energy(
                    name=lig.replace('.mol', ''),
                    mol=stk_lig,
                    ey_file=ey_file
                )
            # Read energy.
            # kJ/mol.
            E_extracted = read_gfnx2xtb_eyfile(ey_file)

            # Calculate energy of optimised ligand.
            # Load in lowest energy conformer.
            opt_mol = stk.BuildingBlock.init_from_file(
                filename_
            )
            if not os.path.exists(opt_lig_ey):
                calculate_energy(
                    name=opt_lig_n,
                    mol=opt_mol,
                    ey_file=opt_lig_ey
                )
            # Read energy.
            # kJ/mol.
            E_free = read_gfnx2xtb_eyfile(opt_lig_ey)
            # Add to list the strain energy:
            # (E(extracted) - E(optimised/free))
            lse = E_extracted - E_free
            # kJ/mol.
            strain_energies[lig] = lse

        # Write data.
        with open(output_json, 'w') as f:
            json.dump(strain_energies, f)

    # Get data.
    with open(output_json, 'r') as f:
        strain_energies = json.load(f)

    return strain_energies


def get_furthest_pair_FGs(stk_mol):
    """
    Returns the pair of functional groups that are furthest apart.

    """

    if stk_mol.get_num_functional_groups() == 2:
        return tuple(i for i in stk_mol.get_functional_groups())
    elif stk_mol.get_num_functional_groups() < 2:
        raise ValueError(f'{stk_mol} does not have at least 2 FGs')

    fg_centroids = [
        (fg, stk_mol.get_centroid(atom_ids=fg.get_placer_ids()))
        for fg in stk_mol.get_functional_groups()
    ]

    fg_dists = sorted(
        [
            (i[0], j[0], euclidean(i[1], j[1]))
            for i, j in combinations(fg_centroids, 2)
        ],
        key=lambda x: x[2],
        reverse=True
    )

    return (fg_dists[0][0], fg_dists[0][1])


def calculate_deltaangle_distance(
    org_ligs,
    smiles_keys,
    fg_factory,
    file_prefix=None
):

    """
    Calculate the change of bite angle of each ligand in the cage.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC and metals + AromaticCNN.

    Parameters
    ----------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    fg_factory :
    :class:`iterable` of :class:`stk.FunctionalGroupFactory`
        Functional groups to asign to molecules.
        NN_distance calculator will not work for cages built from FGs
        other than metals + AromaticCNC and metals + AromaticCNN.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns
    -------
    delta_angles : :class:`dict`
        Bite angle in cage - free optimised ligand for each ligand.
        Output is absolute values.

    """

    delta_angles = {}
    # Iterate over ligands.
    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        smiles_key = stk.Smiles().get_key(stk_lig)
        idx = smiles_keys[smiles_key]
        sgt = str(stk_lig.get_num_atoms())
        # Get optimized ligand name that excludes any cage
        # information.
        if file_prefix is None:
            filename_ = f'organic_linker_s{sgt}_{idx}_opt.mol'
        else:
            filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'

        _in_cage = stk.BuildingBlock.init_from_molecule(
            stk_lig,
            functional_groups=fg_factory
        )
        _in_cage = _in_cage.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(_in_cage)
        )

        _free = stk.BuildingBlock.init_from_file(
            filename_,
            functional_groups=fg_factory
        )
        _free = _free.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(_free)
        )
        angle_in_cage = calculate_bite_angle(bb=_in_cage)
        angle_free = calculate_bite_angle(bb=_free)

        delta_angles[lig] = abs(angle_in_cage - angle_free)

    return delta_angles


def calculate_deltann_distance(
    org_ligs,
    smiles_keys,
    fg_factory,
    file_prefix=None
):
    """
    Calculate the change of NN distance of each ligand in the cage.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC and metals + AromaticCNN.

    Parameters
    ----------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    fg_factory :
    :class:`iterable` of :class:`stk.FunctionalGroupFactory`
        Functional groups to asign to molecules.
        NN_distance calculator will not work for cages built from FGs
        other than metals + AromaticCNC and metals + AromaticCNN.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns
    -------
    delta_nns : :class:`dict`
        NN distance in cage - free optimised ligand for each ligand.
        Output is absolute values.

    """

    delta_nns = {}
    # Iterate over ligands.
    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        smiles_key = stk.Smiles().get_key(stk_lig)
        idx = smiles_keys[smiles_key]
        sgt = str(stk_lig.get_num_atoms())
        # Get optimized ligand name that excludes any cage
        # information.
        if file_prefix is None:
            filename_ = f'organic_linker_s{sgt}_{idx}_opt.mol'
        else:
            filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'

        _in_cage = stk.BuildingBlock.init_from_molecule(
            stk_lig,
            functional_groups=fg_factory
        )
        _in_cage = _in_cage.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(_in_cage)
        )

        _free = stk.BuildingBlock.init_from_file(
            filename_,
            functional_groups=fg_factory
        )
        _free = _free.with_functional_groups(
            functional_groups=get_furthest_pair_FGs(_free)
        )

        nn_in_cage = calculate_NN_distance(bb=_in_cage)
        nn_free = calculate_NN_distance(bb=_free)

        delta_nns[lig] = abs(nn_in_cage - nn_free)

    return delta_nns


def get_organic_linkers(cage, metal_atom_nos, file_prefix=None):
    """
    Extract a list of organic linker .Molecules from a cage.

    Parameters
    ----------
    cage : :class:`stk.Molecule`
        Molecule to get the organic linkers from.

    metal_atom_nos : :class:`iterable` of :class:`int`
        The atomic number of metal atoms to remove from structure.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns
    -------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    """

    org_lig = {}

    # Produce a graph from the cage that does not include metals.
    cage_g = nx.Graph()
    atom_ids_in_G = set()
    for atom in cage.get_atoms():
        if atom.get_atomic_number() in metal_atom_nos:
            continue
        cage_g.add_node(atom)
        atom_ids_in_G.add(atom.get_id())

    # Add edges.
    for bond in cage.get_bonds():
        a1id = bond.get_atom1().get_id()
        a2id = bond.get_atom2().get_id()
        if a1id in atom_ids_in_G and a2id in atom_ids_in_G:
            cage_g.add_edge(bond.get_atom1(), bond.get_atom2())

    # Get disconnected subgraphs as molecules.
    # Sort and sort atom ids to ensure molecules are read by RDKIT
    # correctly.
    connected_graphs = [
        sorted(subgraph, key=lambda a: a.get_id())
        for subgraph in sorted(nx.connected_components(cage_g))
    ]
    smiles_keys = {}
    for i, cg in enumerate(connected_graphs):
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = [i.get_id() for i in atoms]
        cage.write(
            'temporary_linker.mol',
            atom_ids=atom_ids
        )
        temporary_linker = stk.BuildingBlock.init_from_file(
            'temporary_linker.mol'
        ).with_canonical_atom_ordering()
        smiles_key = stk.Smiles().get_key(temporary_linker)
        if smiles_key not in smiles_keys:
            smiles_keys[smiles_key] = len(smiles_keys.values())+1
        idx = smiles_keys[smiles_key]
        sgt = str(len(atoms))
        # Write to mol file.
        if file_prefix is None:
            filename_ = f'organic_linker_s{sgt}_{idx}_{i}.mol'
        else:
            filename_ = f'{file_prefix}{sgt}_{idx}_{i}.mol'

        org_lig[filename_] = temporary_linker
        os.system('rm temporary_linker.mol')
        # Rewrite to fix atom ids.
        org_lig[filename_].write(filename_)
        org_lig[filename_] = stk.BuildingBlock.init_from_file(
            filename_
        )

    return org_lig, smiles_keys


def get_lowest_energy_conformers(
    org_ligs,
    smiles_keys,
    file_prefix=None,
    gfn_exec=None,
    conformer_function=None,
    conformer_settings=None,
):
    """
    Determine the lowest energy conformer of cage organic linkers.

    Will do multiple if there are multiple types.

    Parameters
    ----------
    org_ligs : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    gfn_exec : :class:`str`, optional
        Location of GFN-xTB executable to use.

    conformer_function : :class:`function`, optional
        Define the function used to rank and find the lowest energy
        conformer.

    """

    if conformer_function is None:
        conformer_function = get_lowest_energy_conformer
    if conformer_settings is None:
        conformer_settings = None

    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        smiles_key = stk.Smiles().get_key(stk_lig)
        idx = smiles_keys[smiles_key]
        sgt = str(stk_lig.get_num_atoms())
        # Get optimized ligand name that excludes any cage information.
        if file_prefix is None:
            filename_ = f'organic_linker_s{sgt}_{idx}_opt.mol'
            ligand_name_ = f'organic_linker_s{sgt}_{idx}_opt'
        else:
            filename_ = f'{file_prefix}{sgt}_{idx}_opt.mol'
            ligand_name_ = f'{file_prefix}{sgt}_{idx}_opt'

        if not os.path.exists(filename_):
            if not os.path.exists(f'{ligand_name_}_confs/'):
                os.mkdir(f'{ligand_name_}_confs/')
            low_e_conf = conformer_function(
                name=ligand_name_,
                mol=stk_lig,
                gfn_exec=gfn_exec,
                settings=conformer_settings
            )
            low_e_conf.write(filename_)


class MissingSettingError(Exception):
    ...


def build_conformers(mol, N, ETKDG_version=None):
    """
    Convert stk mol into RDKit mol with N conformers.

    ETKDG_version allows the user to pick their choice of ETKDG params.

    `None` provides the settings used in ligand_combiner and unsymm.

    Other options:
        `v3`:
            New version from DOI: 10.1021/acs.jcim.0c00025
            with improved handling of macrocycles.

    """
    molecule = mol.to_rdkit_mol()
    molecule.RemoveAllConformers()

    if ETKDG_version is None:
        cids = rdkit.EmbedMultipleConfs(
            mol=molecule,
            numConfs=N,
            randomSeed=1000,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
            numThreads=4,
        )

    elif ETKDG_version == 'v3':
        params = rdkit.ETKDGv3()
        params.randomSeed = 1000
        cids = rdkit.EmbedMultipleConfs(
            mol=molecule,
            numConfs=N,
            params=params
        )

    print(f'there are {molecule.GetNumConformers()} conformers')
    return cids, molecule


def get_lowest_energy_conformer(
    name,
    mol,
    gfn_exec=None,
    settings=None,
):
    """
    Get lowest energy conformer of molecule.

    Method:
        1) ETKDG conformer search on molecule
        2) xTB `normal` optimisation of each conformer
        3) xTB `opt_level` optimisation of lowest energy conformer
        4) save file

    """

    if settings is None:
        settings = {
            'conf_opt_level': 'normal',
            'final_opt_level': 'extreme',
            'charge': 0,
            'no_unpaired_e': 0,
            'max_runs': 1,
            'calc_hessian': False,
            'solvent': None,
            'N': 100
        }

    # Check for missing settings.
    req_settings = [
        'N', 'final_opt_level', 'charge', 'no_unpaired_e', 'max_runs',
        'calc_hessian', 'solvent', 'conf_opt_level'
    ]
    for i in req_settings:
        if i not in settings:
            raise MissingSettingError(
                f'Settings missing {i}. Has {settings.keys()}.'
            )

    # Run ETKDG on molecule.
    print(f'....running ETKDG on {name}')
    cids, confs = build_conformers(mol, N=settings['N'])

    # Optimize all conformers at normal level with xTB.
    low_e_conf_id = -100
    low_e = 10E20
    for cid in cids:
        name_ = f'{name}_confs/c_{cid}'
        ey_file = f'{name}_confs/c_{cid}_eyout'

        mol = update_from_rdkit_conf(
            mol,
            confs,
            conf_id=cid
        )
        mol.write(f'{name}_confs/c_{cid}.mol')

        # Optimize.
        opt_mol = optimize_conformer(
            name=name_,
            mol=mol,
            gfn_exec=gfn_exec,
            opt_level=settings['conf_opt_level'],
            charge=settings['charge'],
            no_unpaired_e=settings['no_unpaired_e'],
            max_runs=settings['max_runs'],
            calc_hessian=settings['calc_hessian'],
            solvent=settings['solvent']
        )
        opt_mol.write(f'{name}_confs/c_{cid}_opt.mol')

        # Get energy.
        calculate_energy(
            name=name_,
            mol=opt_mol,
            gfn_exec=gfn_exec,
            ey_file=ey_file,
            charge=settings['charge'],
            no_unpaired_e=settings['no_unpaired_e'],
            solvent=settings['solvent']
        )
        ey = read_gfnx2xtb_eyfile(ey_file)
        if ey < low_e:
            print(
                'lowest energy conformer updated with energy: '
                f'{ey}, id: {cid}'
            )
            low_e_conf_id = cid
            low_e = ey

    # Get lowest energy conformer.
    low_e_conf = stk.BuildingBlock.init_from_file(
        f'{name}_confs/c_{low_e_conf_id}_opt.mol'
    )
    low_e_conf.write(f'{name}_confs/low_e_unopt.mol')

    # Optimize lowest energy conformer at opt_level.
    low_e_conf = optimize_conformer(
        name=name+'low_e_opt',
        mol=low_e_conf,
        gfn_exec=gfn_exec,
        opt_level=settings['final_opt_level'],
        charge=settings['charge'],
        no_unpaired_e=settings['no_unpaired_e'],
        max_runs=settings['max_runs'],
        calc_hessian=settings['calc_hessian'],
        solvent=settings['solvent']
    )
    low_e_conf.write(f'{name}_confs/low_e_opt.mol')

    # Return molecule.
    return low_e_conf


def optimize_conformer(
    name,
    mol,
    gfn_exec=None,
    opt_level='extreme',
    charge=0,
    no_unpaired_e=0,
    max_runs=1,
    calc_hessian=False,
    solvent=None
):
    """
    Run simple GFN-xTB optimisation of molecule.

    """

    if gfn_exec is None:
        gfn_exec = '/home/atarzia/software/xtb-190806/bin/xtb'

    print(f'....optimizing {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    xtb_opt = stko.XTB(
        xtb_path=gfn_exec,
        output_dir=f'{name}_opt',
        gfn_version=2,
        num_cores=6,
        opt_level=opt_level,
        charge=charge,
        num_unpaired_electrons=no_unpaired_e,
        max_runs=max_runs,
        calculate_hessian=calc_hessian,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )

    return xtb_opt.optimize(mol=mol)
