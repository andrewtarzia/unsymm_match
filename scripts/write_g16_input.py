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
import stk


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


def write_defn_section(
    name,
    runtype,
    restart=False,
    solvent=None,
):

    if solvent is None:
        solvent_s = ''
    else:
        solvent_s = solvent

    new_name = name.split('.')[0]

    if restart:
        geom_ = 'geom=(checkpoint)'
    else:
        geom_ = ''

    string = (
        '#P PBE1PBE/GenECP '
        f'{runtype} '
        f'{geom_} '
        'SCF=(YQC,MaxCycle=900) '
        'int=(Grid=Superfinegrid) '
        'EmpiricalDispersion=GD3BJ '
        f'{solvent_s}\n\n'
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


def write_run_file(name, directory, infile, np):

    runfile = f'{directory}/{name}.sh'

    runlines = f"g16 < {infile} > {infile.replace('.gau', '.log')}\n"

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


def write_opt_input_file(
    infile,
    struct,
    np,
    directory,
    runtype,
    org_basis,
    m_basis,
    solvent,
):
    base_name = '_'+infile.replace('.in', '')

    string = write_top_section(
        name=base_name,
        np=np,
        restart=False
    )
    string += write_defn_section(
        name=base_name,
        runtype=runtype,
        restart=False,
        solvent=solvent,
    )
    string += write_molecule_section(struct, restart=False)
    string += write_ECP_section(org_basis=org_basis, m_basis=m_basis)

    string2 = write_top_section(
        name=base_name,
        np=np,
        restart=True
    )
    string2 += write_defn_section(
        name=base_name,
        runtype=runtype,
        restart=True,
        solvent=solvent,
    )
    string2 += write_molecule_section(struct, restart=True)
    string2 += write_ECP_section(org_basis=org_basis, m_basis=m_basis)

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
        runtype='SP',
        restart=True,
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
        runtype='freq',
        restart=True,
        freq=True,
    )
    string += write_molecule_section(struct, restart=True)
    string += write_ECP_section(org_basis='Def2SVP', m_basis='SDD')

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def main():
    num_proc = 32
    optc_extension = '_optc.mol'
    isomers = ['A', 'B', 'C', 'D']
    selected_ligands = [
        '5D1', '4D2', '5D3', '3D1',
        '4B3', '4B1', '5B4',
        '5A3', '5A1',
        '4C1', '4C3',
        # '5C2', '4C2', '3C2',
    ]

    # full_data = pd.read_csv('../all_cage_results.txt')
    # df_to_run = full_data[full_data['lig'].isin(selected_ligands)]

    for lig in selected_ligands:
        # row = df_to_run[df_to_run['lig'] == lig]
        for isomer in isomers:
            # param = (
            #     float(row[f'energy_{isomer}']),
            #     float(row[f'plane_dev_{isomer}']),
            #     float(row[f'sqpl_op_{isomer}']),
            # )
            struct_file = f'../{lig}_{isomer}{optc_extension}'
            prefix = f'{lig.lower()}_{isomer.lower()}'
            struct = stk.BuildingBlock.init_from_file(struct_file)
            opt_directory = f'opt_{prefix}'
            if not exists(opt_directory):
                mkdir(opt_directory)

            step1_infile = f'o1_{prefix}.gau'
            write_opt_input_file(
                infile=step1_infile,
                struct=struct,
                np=num_proc,
                directory=opt_directory,
                runtype='opt=(calcfc,maxstep=4)',
                org_basis='Def2SVP',
                m_basis='LANL2DZ',
                solvent=None,
            )
            write_run_file(
                name=f'o1_{prefix}',
                directory=opt_directory,
                infile=step1_infile,
                np=num_proc
            )

            step2_infile = f'o2_{prefix}.gau'
            write_opt_input_file(
                infile=step2_infile,
                struct=struct,
                np=num_proc,
                directory=opt_directory,
                runtype='opt=(maxstep=4)',
                org_basis='Def2SVP',
                m_basis='SDD',
                solvent=None,
            )
            write_run_file(
                name=f'o2_{prefix}',
                directory=opt_directory,
                infile=step2_infile,
                np=num_proc
            )

            step3_infile = f'o3_{prefix}.gau'
            write_opt_input_file(
                infile=step3_infile,
                struct=struct,
                np=num_proc,
                directory=opt_directory,
                runtype='opt=(maxstep=4)',
                org_basis='Def2SVP',
                m_basis='SDD',
                solvent=r'SCRF=(PCM,Solvent=DiMethylSulfoxide)',
            )
            write_run_file(
                name=f'o3_{prefix}',
                directory=opt_directory,
                infile=step3_infile,
                np=num_proc
            )


if __name__ == '__main__':
    main()
