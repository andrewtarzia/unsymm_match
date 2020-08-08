#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write input for ORCA 4.2.1 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

import json
from os.path import exists
from os import mkdir
import subprocess as sp
import glob
import stk


def load_in_structure(file):
    struct = stk.BuildingBlock.init_from_file(file)
    pos_mat = struct.get_position_matrix()

    mol_list = []
    for i, atom in enumerate(struct.get_atoms()):
        mol_list.append([
            atom.__class__.__name__,
            (pos_mat[i][0], pos_mat[i][1], pos_mat[i][2])
        ])

    return mol_list


def write_molecule_section(mol_list):

    charge = 4
    multiplicity = 1

    string = f'* xyz {charge} {multiplicity}\n'
    for atom in mol_list:
        string += (
            f'{atom[0]} {round(atom[1][0], 4)} '
            f'{round(atom[1][1], 4)} {round(atom[1][2], 4)}\n'
        )
    string += '*\n'

    return string


def write_spe_input_file(infile, mol_list, grid, np, directory):
    comment_line = (
        '# Test a simple DFT calculation with D3, solvent (DMSO, CPCM)'
        ', ECPs. Using RIJCOSX because PBE0 is a hybrid functional.\n'
    )
    top_line = (
        f'! DFT SP RKS PBE0 RIJCOSX def2-SVP D3BJ {grid} NOFINALGRID '
        'TightSCF CPCM(dmso) printbasis\n\n'
    )
    base_name = infile.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 1000\nend\n\n'
    )

    procs_section = (
        f'%pal\n   nprocs {np}\nend\n\n'
    )

    mol_section = write_molecule_section(mol_list)

    string = comment_line
    string += top_line
    string += base_line
    string += scf_section
    string += procs_section
    string += mol_section

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def write_run_file(name, directory, infiles, np):

    orca_bin_dir = '/apps/orca/4.2.1/bin/orca'

    runfile = f'{directory}/{name}.sh'

    runlines = ''.join([
        f"{orca_bin_dir} {infile} > {infile.replace('.in', '.out')}\n"
        for infile in infiles
    ])

    string = (
        f'#PBS -N {name}\n'
        '#PBS -l walltime=72:00:00\n'
        f'#PBS -l select=1:ncpus={np}:mem=124gb\n\n'
        'module load orca/4.2.1\n\n'
        'cd $PBS_O_WORKDIR\n\n'
        f'{runlines}'
    )

    with open(runfile, 'w') as f:
        f.write(string)


def write_opt_input_file(infile, mol_list, grid):
    raise NotImplementedError()
    comment_line = (
        '# Test a simple DFT calculation with D3, solvent (DMSO, CPCM)'
        ', ECPs. Using RIJCOSX because PBE0 is a hybrid functional.\n'
    )
    top_line = (
        f'! DFT SP RKS PBE0 RIJCOSX def2-SVP D3BJ {grid} NOFINALGRID '
        'TightSCF CPCM(dmso) printbasis\n\n'
    )
    base_name = infile.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 1000\nend\n'
    )

    mol_section = write_molecule_section(mol_list)

    string = (
        comment_line + top_line + base_line + scf_section + mol_section
    )

    with open(infile, 'w') as f:
        f.write(string)


def get_final_energy(outfile):

    with open(outfile, 'r') as f:
        for line in f.readlines():
            if 'FINAL SINGLE POINT ENERGY' in line:
                energy = float(line.rstrip().split(' ')[-1])
                print(line, energy)
                break

    return energy


def run_spe(calc_name, mol_list, grid):

    infile = f'{calc_name}.in'
    outfile = infile.replace('.in', '.out')
    write_spe_input_file(infile, mol_list, grid)

    orca_dir = '/home/atarzia/software/orca/orca'
    cmd = f'{orca_dir} {infile}'

    with open(outfile, 'w') as f:
        # Note that sp.call will hold the program until completion
        # of the calculation.
        sp.call(
            cmd,
            stdin=sp.PIPE,
            stdout=f,
            stderr=sp.PIPE,
            # Shell is required to run complex arguments.
            shell=True
        )

    energy = get_final_energy(outfile)

    return energy


def output_energies(energy_file, energy_dict):

    with open(energy_file, 'w') as f:
        json.dump(energy_dict, f)


def main():

    list_of_mols = glob.glob('*.mol')
    grid_lists = [
        'Grid0', 'Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5', 'Grid6',
        'Grid7',
    ]
    num_proc = 32

    for moll in list_of_mols:
        prefix = moll.replace('.mol', '')

        mol_list = load_in_structure(moll)
        infiles = []
        for grid in grid_lists:
            calc_name = f'{prefix}_spe_{grid}'
            directory = f'dir_{prefix}'
            if not exists(directory):
                mkdir(directory)
            print(f'> writing {calc_name}.....')
            infile = f'{calc_name}.in'
            write_spe_input_file(
                infile=infile,
                mol_list=mol_list,
                grid=grid,
                np=num_proc,
                directory=directory,
            )
            infiles.append(infile)
            continue
            write_opt_input_file(infile, mol_list, grid)
        write_run_file(
            name=prefix,
            directory=directory,
            infiles=infiles,
            np=num_proc
        )


if __name__ == '__main__':
    main()
