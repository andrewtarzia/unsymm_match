#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write input for Gaussian 16 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

from os.path import exists
from os import mkdir
import glob
import stk


def load_in_structure(file):

    return stk.BuildingBlock.init_from_file(file)


def write_top_section(name, np, restart=False, freq=False):

    new_name = name.split('.')[0]

    if freq:
        string = (
            f'%RWF={new_name}.rwf\n'
            f'%NoSave\n'
            f'%oldchk={new_name}.chk\n'
            f'%chk={new_name}_r.chk\n'
            f'%mem=60000mb\n'
            f'%nprocshared={np}\n\n'
        )
    elif restart:
        string = (
            f'%oldchk={new_name}.chk\n'
            f'%chk={new_name}_r.chk\n'
            f'%mem=60000mb\n'
            f'%nprocshared={np}\n\n'
        )
    else:
        string = (
            f'%chk={new_name}.chk\n'
            f'%mem=60000mb\n'
            f'%nprocshared={np}\n\n'
        )

    return string


def write_defn_section(name, freq=False, restart=False, spe=False):

    new_name = name.split('.')[0]

    if freq:
        string = (
            '#P PBE1PBE/GenECP '
            'freq '
            'geom=(checkpoint) '
            'SCF=(YQC,MaxCycle=900) '
            'int=(Grid=Ultrafine) '
            'EmpiricalDispersion=GD3 '
            'SCRF=(PCM,Solvent=DiMethylSulfoxide)\n\n'
            f'{new_name}\n\n'
        )
    elif spe:
        string = (
            '#P PBE1PBE/GenECP '
            'SP '
            'geom=(checkpoint) '
            'SCF=(YQC,MaxCycle=900) '
            'int=(Grid=Ultrafine) '
            'EmpiricalDispersion=GD3 '
            'SCRF=(PCM,Solvent=DiMethylSulfoxide)\n\n'
            f'{new_name}\n\n'
        )
    elif restart:
        string = (
            '#P PBE1PBE/GenECP '
            'opt=(Tight,calcfc) '
            'geom=(checkpoint) '
            'SCF=(YQC,MaxCycle=900) '
            'int=(Grid=Ultrafine) '
            'EmpiricalDispersion=GD3 '
            'SCRF=(PCM,Solvent=DiMethylSulfoxide)\n\n'
            f'{new_name}\n\n'
        )
    else:
        string = (
            '#P PBE1PBE/GenECP '
            'opt=(Tight,calcfc) '
            'SCF=(YQC,MaxCycle=900) '
            'int=(Grid=Ultrafine) '
            'EmpiricalDispersion=GD3 '
            'SCRF=(PCM,Solvent=DiMethylSulfoxide)\n\n'
            f'{new_name}\n\n'
        )

    return string


def write_ECP_section(org_basis, m_basis):

    string = (
        f'C H N 0\n{org_basis}\n****\n'
        f'Pd 0\n{m_basis}\n****\n'
        '\n'
        f'Pd 0\n{m_basis}\n\n\n\n'
    )

    return string


def write_molecule_section(struct, restart=False):

    charge = 4
    multiplicity = 1

    if restart:
        string = f'{charge} {multiplicity}\n\n'
    else:
        pos_mat = struct.get_position_matrix()
        string = f'{charge} {multiplicity}\n'
        for atom in struct.get_atoms():
            E = atom.__class__.__name__
            x, y, z = pos_mat[atom.get_id()]
            string += (
                f'{E} {round(x, 4)} {round(y, 4)} {round(z, 4)}\n'
            )
        string += '\n'

    return string


def write_run_file(name, directory, infiles, np):

    runfile = f'{directory}/{name}.sh'

    runlines = ''.join([
        f"g16 < {infile} > {infile.replace('.gau', '.log')}\n"
        for infile in infiles
    ])

    string = (
        f'#PBS -N _{name}\n'
        f'#PBS -l select=1:ncpus={np}:mem=124gb\n'
        '#PBS -l walltime=72:00:00\n'
        'module load gaussian/g16-c01-avx\n\n'
        'cd $PBS_O_WORKDIR\n\n'
        f'{runlines}'
    )

    with open(runfile, 'w') as f:
        f.write(string)


def write_opt_input_file(infile, struct, np, directory):
    base_name = '_'+infile.replace('.in', '')

    string = write_top_section(
        name=base_name,
        np=np,
        restart=False
    )
    string += write_defn_section(
        name=base_name,
        restart=False,
        freq=False,
        spe=False,
    )
    string += write_molecule_section(struct, restart=False)
    string += write_ECP_section(org_basis='Def2SVP', m_basis='SDD')

    string2 = write_top_section(
        name=base_name,
        np=np,
        restart=True
    )
    string2 += write_defn_section(
        name=base_name,
        restart=True,
        freq=False,
        spe=False,
    )
    string2 += write_molecule_section(struct, restart=True)
    string2 += write_ECP_section(org_basis='Def2SVP', m_basis='SDD')

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)
    res_infile = infile.replace('.gau', 'r.gau')
    with open(f'{directory}/{res_infile}', 'w') as f:
        f.write(string2)


def write_final_spe_input_file(infile, struct, np, directory):
    base_name = '_'+infile.replace('.in', '')

    string = write_top_section(
        name=base_name,
        np=np,
        restart=True,
    )
    string += write_defn_section(
        name=base_name,
        restart=False,
        freq=False,
        spe=True,
    )
    string += write_molecule_section(struct, restart=True)
    string += write_ECP_section(org_basis='Def2TZVP', m_basis='SDD')

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def write_freq_input_file(infile, struct, np, directory):
    base_name = '_'+infile.replace('.in', '')

    string = write_top_section(
        name=base_name,
        np=np,
        restart=False,
        freq=True,
    )
    string += write_defn_section(
        name=base_name,
        restart=False,
        freq=True,
        spe=False,
    )
    string += write_molecule_section(struct, restart=True)
    string += write_ECP_section(org_basis='Def2SVP', m_basis='SDD')

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def main():

    list_of_mols = glob.glob('*.mol')
    num_proc = 32

    for moll in list_of_mols:
        prefix = moll.replace('.mol', '')

        struct = load_in_structure(moll)
        opt_infiles = []
        fin_spe_infiles = []
        freq_infiles = []
        # for grid in grid_lists:
        calc_name = f'{prefix}'
        opt_directory = f'opt_{prefix}'
        if not exists(opt_directory):
            mkdir(opt_directory)
        print(f'> writing {calc_name}.....')
        o_infile = f'o_{calc_name}.gau'
        s_infile = f's_{calc_name}.gau'
        f_infile = f'f_{calc_name}.gau'
        write_opt_input_file(
            infile=o_infile,
            struct=struct,
            np=num_proc,
            directory=opt_directory,
        )
        opt_infiles.append(o_infile)
        write_final_spe_input_file(
            infile=s_infile,
            struct=struct,
            np=num_proc,
            directory=opt_directory,
        )
        fin_spe_infiles.append(s_infile)
        write_freq_input_file(
            infile=f_infile,
            struct=struct,
            np=num_proc,
            directory=opt_directory,
        )
        freq_infiles.append(f_infile)
        write_run_file(
            name=f'o_{prefix}',
            directory=opt_directory,
            infiles=opt_infiles,
            np=num_proc
        )
        write_run_file(
            name=f's_{prefix}',
            directory=opt_directory,
            infiles=fin_spe_infiles,
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
