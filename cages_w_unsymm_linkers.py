#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all isomers of a MOC with nonsymmetric linkers.

Author: Andrew Tarzia

Date Created: 22 Jun 2019

"""

import sys
import logging
import stk


def edge_dirs():
    """
    Returns list of edge directions for all isomers.

    Based on stk.M2L4_Lantern() with two building blocks and includes
    the heteroleptic cages.

    Returns
    -------
    :class:`list`
        List of building block placements to define cage.

    """
    return [
        [1, 1, 1, 1],
        [1, -1, 1, 1],
        [1, -1, 1, -1],
        [1, 1, -1, -1],
    ]


def main():
    if (not len(sys.argv) == 2):
        print("""
    Usage: cages_w_unsymm_linkers.py linker
        linker (str) - .mol file containing cage linker
        """)
        sys.exit()
    else:
        linker = sys.argv[1]

    linker_su = stk.StructUnit2(linker, ['pyridine_N_metal'])

    bb_dirs = edge_dirs()

    # From input structure, determine precursors and calculate their
    # energies for use later.
    # Define metal StructUnit and dir to precursor file.
    metal_su = stk.MetalStructUnit(
        element='Pd',
        coord_no=4,
        geometry='sqpl',
        charge=2
    )

    # Iterate through cages and calculate formation energies.
    for bb in bb_dirs:
        logging.info(
            f'Doing isomer with ordering {bb}.'
        )
        # Build cage.
        topology =stk.M2L4_Lantern(
            vertex_maps={
                # vertices to put metals
                metal_su: [0, 1]
            },
            edge_maps={
                # edges to put linkers
                linker_su: [0, 1, 2, 3]
            },
            vertex_aligners=[0, 0, 0, 0],
            edge_directions=bb
        )
        cage = stk.MOC(
            building_blocks=[
                metal_su,
                linker_su
            ],
            topology=topology
        )

        # Output to files.
        file_prefix = f"MOC_{''.join(str(i) for i in bb)}"
        cage.write(f'{file_prefix}.mol')
        cage.write(f'{file_prefix}.pdb')
        cage.dump(f'{file_prefix}.json')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
