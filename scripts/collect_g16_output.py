#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to collect output from Gaussian 16 DFT SPE calculations.

Author: Andrew Tarzia

Date Created: 5 Aug 2020

"""

from os.path import join
from os import walk
import json
import glob


def get_g16_energy(filename, front_splitter, back_splitter):

    with open(filename, 'r') as f:
        data = f.read()
    energy = float(
        data.split(front_splitter)[2].split(back_splitter)[0]
    )

    return energy  # a.u.


def main():

    collated_energies = {}
    for dir in next(walk('.'))[1]:
        print(dir)
    with open('collated_energies.json', 'w') as f:
        json.dump(collated_energies, f)


if __name__ == '__main__':
    main()
