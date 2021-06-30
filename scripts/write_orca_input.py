#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to write input for Orca 4.2 DFT calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

from os.path import exists
from os import mkdir
import stk


def write_top_section(basename, method, solvent):

    if solvent is None:
        solvent_s = ''
    else:
        solvent_s = solvent

    if method == 'PBE0':
        string = (
            f'! DFT SP RKS PBE0 def2-SVP D4 Grid6 NOFINALGRID '
            f'SlowConv TightSCF printbasis {solvent_s}'
            '\n\n'
            f'%base "{basename}"\n'
            '%maxcore 3000\n'
            '%scf\n   MaxIter 2000\nend\n'
        )
    elif method == 'B97-3c':
        string = (
            f'! DFT SP B97-3c TightSCF printbasis {solvent_s} Grid6 '
            'NOFINALGRID SlowConv '
            '\n\n'
            f'%base "{basename}"\n'
            '%maxcore 3000\n'
            '%scf\n   MaxIter 2000\nend\n'
        )

    return string


def write_proc_section(np):

    return f'%pal\n   nprocs {np}\nend\n\n'


def write_molecule_section(filename):

    charge = 4
    multi = 1

    return f'* xyzfile {charge} {multi} {filename}\n'


def write_input_file(
    prefix,
    np,
    directory,
    method,
    solvent,
):

    basename = f'_s_{prefix}_{method}'
    infile = f's_{prefix}_{method}.in'

    string = write_top_section(
        basename=basename,
        method=method,
        solvent=solvent,
    )
    string += write_proc_section(np)
    string += write_molecule_section(f's_{prefix}.xyz')

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

    methods = {'pbe': 'PBE0', 'b973c': 'B97-3c'}

    opt_directory = 's_orca'
    if not exists(opt_directory):
        mkdir(opt_directory)
    for lig in selected_ligands:
        for isomer in isomers:
            struct_file = f'../{lig}_{isomer}{optc_extension}'
            prefix = f'{lig.lower()}_{isomer.lower()}'
            struct = stk.BuildingBlock.init_from_file(struct_file)
            struct.write(f'{opt_directory}/s_{prefix}.xyz')

            for meth in methods:
                write_input_file(
                    prefix=prefix,
                    np=num_proc,
                    directory=opt_directory,
                    method=methods[meth],
                    solvent=None,
                )


if __name__ == '__main__':
    main()
