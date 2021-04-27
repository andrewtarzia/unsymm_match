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
import glob
import stk


def load_in_structure(file):
    struct = stk.BuildingBlock.init_from_file(file)
    # pos_mat = struct.get_position_matrix()

    # mol_list = []
    # for i, atom in enumerate(struct.get_atoms()):
    #     mol_list.append([
    #         atom.__class__.__name__,
    #         (pos_mat[i][0], pos_mat[i][1], pos_mat[i][2])
    #     ])

    return struct


def write_molecule_section(directory, base_name, struct):

    charge = 4
    multiplicity = 1
    xyzfile = f'{base_name}_init.xyz'

    string = f'* xyzfile {charge} {multiplicity} {xyzfile}\n'
    # for atom in mol_list:
    #     string += (
    #         f'{atom[0]} {round(atom[1][0], 4)} '
    #         f'{round(atom[1][1], 4)} {round(atom[1][2], 4)}\n'
    #     )
    # string += '*\n'
    struct.write(f'{directory}/{xyzfile}')

    return string


def write_spe_input_file(infile, struct, grid, np, directory):
    comment_line = (
        '# DFT single point energy with D4, solvent (DMSO, CPCM), '
        'def2-TZVPP, def2-ECPs. Grid is 6 and slow convergence is on.'
        '.\n'
    )
    top_line = (
        f'! DFT SP RKS PBE0 def2-TZVPP def2/J D4 '
        f'Grid{grid} NOFINALGRID SlowConv '
        'TightSCF CPCM(dmso) printbasis\n\n'
    )
    base_name = '_'+infile.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 2000\nend\n\n'
    )

    procs_section = (
        f'%pal\n   nprocs {np}\nend\n\n'
    )

    mol_section = write_molecule_section(directory, base_name, struct)

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


def write_opt_input_file(infile, struct, grid, np, directory):
    comment_line = (
        '# DFT optimisation with D3, solvent (DMSO, CPCM)'
        ', def2-SVP, def2-ECPs. Grid is 6 and slow convergence is on.'
        'Using RIJCOSX with GridX6 because PBE0 is a hybrid functional'
        '.\n'
    )
    top_line = (
        f'! DFT OPT RKS PBE0 RIJCOSX GridX{grid} def2-SVP def2/J D3BJ '
        f'Grid{grid} NOFINALGRID SlowConv '
        'TightSCF CPCM(dmso) XYZFILE PDBFILE printbasis\n\n'
    )
    base_name = '_'+infile.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 2000\nend\n\n'
    )

    geom_section = (
        '%geom\n   MaxIter 1000\nend\n\n'
    )

    procs_section = (
        f'%pal\n   nprocs {np}\nend\n\n'
    )

    mol_section = write_molecule_section(directory, base_name, struct)

    string = comment_line
    string += top_line
    string += base_line
    string += scf_section
    string += geom_section
    string += procs_section
    string += mol_section

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def write_final_opt_input_file(infile, struct, grid, np, directory):
    comment_line = (
        '# DFT optimisation with D3, solvent (DMSO, CPCM), def2-TZVPP,'
        ' def2-ECPs. Grid is 6 and slow convergence is on.'
        '.\n'
    )
    top_line = (
        f'! DFT OPT RKS PBE0 def2-TZVPP def2/J D3BJ '
        f'Grid{grid} NOFINALGRID SlowConv '
        'TightSCF CPCM(dmso) XYZFILE PDBFILE printbasis\n\n'
    )
    base_name = '_'+infile.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 2000\nend\n\n'
    )

    geom_section = (
        '%geom\n   MaxIter 1000\nend\n\n'
    )

    procs_section = (
        f'%pal\n   nprocs {np}\nend\n\n'
    )

    mol_section = write_molecule_section(directory, base_name, struct)

    string = comment_line
    string += top_line
    string += base_line
    string += scf_section
    string += geom_section
    string += procs_section
    string += mol_section

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def write_freq_input_file(infile, struct, grid, np, directory):
    comment_line = (
        '# DFT frequency calculation with D3, solvent (DMSO, CPCM), '
        'def2-TZVPP, def2-ECPs. Grid is 6 and slow convergence is on.'
        '.\n'
    )
    top_line = (
        f'! DFT NumFreq RKS PBE0 def2-TZVPP def2/J D3BJ '
        f'Grid{grid} NOFINALGRID SlowConv '
        'TightSCF CPCM(dmso) XYZFILE PDBFILE printbasis\n\n'
    )
    base_name = '_'+infile.replace('.in', '')
    base_line = f'%base "{base_name}" \n\n'

    scf_section = (
        '%scf\n   MaxIter 2000\nend\n\n'
    )

    procs_section = (
        f'%pal\n   nprocs {np}\nend\n\n'
    )

    mol_section = write_molecule_section(directory, base_name, struct)

    string = comment_line
    string += top_line
    string += base_line
    string += scf_section
    string += procs_section
    string += mol_section

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def get_spe_energy(dir, basename):

    out_file = f'{dir}/{basename}.out'

    with open(out_file, 'r') as f:
        lines = f.readlines()

    targ_str = 'FINAL SINGLE POINT ENERGY'
    for line in lines:
        if targ_str in line:
            energy = float(line.rstrip().split(' ')[-1])
            break

    return energy * 2625.5


def get_opt_energy(dir, basename):

    out_file = f'{dir}/{basename}_fin.out'

    if not exists(out_file):
        raise FileNotFoundError(f'Copy final out file to {out_file}')

    with open(out_file, 'r') as f:
        lines = f.readlines()

    targ_str = 'FINAL SINGLE POINT ENERGY'
    for line in lines:
        if targ_str in line:
            energy = float(line.rstrip().split(' ')[-1])
            print(energy)

    print(energy)
    input()
    return energy * 2625.5


def get_dft_spe_energies():

    files = glob.glob('*_opt.mol')
    spe_directory = 'spe_calcs'

    energies = {}
    for file in sorted(files):
        print('dft_energy:', file)
        name = file.replace('_opt.mol', '')
        o_output = f'{name}_spe.out'
        basename = f'{name}_spe'
        if exists(f'{spe_directory}/{o_output}'):
            energies[name] = get_spe_energy(spe_directory, basename)
        else:
            energies[name] = -5

    return energies


def get_dft_opt_energies():

    files = glob.glob('*_opt.mol')
    opt_directory = 'opt_calcs'

    energies = {}
    for file in sorted(files):
        print('dft_energy:', file)
        name = file.replace('_opt.mol', '')
        o_output = f'{name}_opt.out'
        basename = f'{name}_opt'
        if exists(f'{opt_directory}/{o_output}'):
            energies[name] = get_opt_energy(opt_directory, basename)
        else:
            energies[name] = -5

    return energies


def output_energies(energy_file, energy_dict):

    with open(energy_file, 'w') as f:
        json.dump(energy_dict, f)


def main():

    list_of_mols = glob.glob('*.mol')
    # grid_lists = [0, 1, 2, 3, 4, 5, 6, 7]
    grid_chosen = 6
    num_proc = 32

    for moll in list_of_mols:
        prefix = moll.replace('.mol', '')

        struct = load_in_structure(moll)
        infiles = []
        opt_infiles = []
        fin_opt_infiles = []
        freq_infiles = []
        # for grid in grid_lists:
        calc_name = f'{prefix}'
        spe_directory = f'spe_{prefix}'
        opt_directory = f'opt_{prefix}'
        if not exists(spe_directory):
            mkdir(spe_directory)
        if not exists(opt_directory):
            mkdir(opt_directory)
        print(f'> writing {calc_name}.....')
        infile = f's_{calc_name}.in'
        o_infile = f'o_{calc_name}.in'
        fo_infile = f'fo_{calc_name}.in'
        f_infile = f'f_{calc_name}.in'
        write_spe_input_file(
            infile=infile,
            struct=struct,
            grid=grid_chosen,
            np=num_proc,
            directory=spe_directory,
        )
        infiles.append(infile)
        # if grid == grid_chosen:
        write_opt_input_file(
            infile=o_infile,
            struct=struct,
            grid=grid_chosen,
            np=num_proc,
            directory=opt_directory,
        )
        opt_infiles.append(o_infile)
        write_final_opt_input_file(
            infile=fo_infile,
            struct=struct,
            grid=grid_chosen,
            np=num_proc,
            directory=opt_directory,
        )
        fin_opt_infiles.append(fo_infile)
        write_freq_input_file(
            infile=f_infile,
            struct=struct,
            grid=grid_chosen,
            np=num_proc,
            directory=opt_directory,
        )
        freq_infiles.append(f_infile)
        write_run_file(
            name=prefix,
            directory=spe_directory,
            infiles=infiles,
            np=num_proc
        )
        write_run_file(
            name=f'o_{prefix}',
            directory=opt_directory,
            infiles=opt_infiles,
            np=num_proc
        )
        write_run_file(
            name=f'fo_{prefix}',
            directory=opt_directory,
            infiles=fin_opt_infiles,
            np=num_proc
        )
        write_run_file(
            name=f'f_{prefix}',
            directory=opt_directory,
            infiles=freq_infiles,
            np=num_proc
        )


if __name__ == '__main__':
    main()
