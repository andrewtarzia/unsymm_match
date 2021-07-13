unsymm_match
============

:author: Andrew Tarzia

Scripts for computational screening ligands for unsymmetrical Pd_2L_4 synthesis.

This work was done while I was a PostDoc in the Jelfs group at Imperial College London <https://github.com/JelfsMaterialsGroup>; <http://www.jelfs-group.org/>.

This work resulted in a publication at chemrxiv <10.26434/chemrxiv.14604294>.

Please contact me with any questions (<atarzia@ic.ac.uk>) or submit an issue!

chemrxiv version: 

.. image:: https://zenodo.org/badge/190806715.svg
   :target: https://zenodo.org/badge/latestdoi/190806715
   
revised version:

.. image:: https://zenodo.org/badge/190806715.svg
   :target: https://zenodo.org/badge/latestdoi/190806715

Installation
------------

To get ``unsymm_match``, you can clone this repo.

Dependancies:

* matplotlib
* pandas
* pymatgen 2019.5.8
* networkx 2.5.1
* mendeleev 0.6.0
* Originally, my package `atools` was used. However, I removed it to avoid a dependency. If a function is missing, the function should be in this release <https://github.com/andrewtarzia/atools/releases/tag/v0.1.unsymm>. Latest version of atools should work, however, the API is ensured to be the same in that release.

`stk`:

    $ pip install stk

`rdkit`, which is a dependency of stk (version 2020.03.1.0 was used in):

    $ conda install -c rdkit rdkit=2020

`stko`:

    $ pip install stko

`xtb`, I used the Version 6.2 RC2 (SAW190805); (190806 binary) throughout this work. However, for future use, a new binary or conda installation of xTB should work. See details at <https://xtb-docs.readthedocs.io/en/latest/contents.html>.

Notes
-----

`screening_process.py` runs the full screening from ligand assembly to cage assembly to analysis and ranking. This script includes all calculations and so should be run in an environment that will not kill it (i.e. using `nohup` on linux).

`ligand_building.py` contains the definition of the ligand components as SMILES strings for editting.

`params` contains one parameter, N, the number of conformers to use in flexibility analysis, which was not used in the manuscript.

`cage_building.py` contains the optimisation protocol that I found to be most successful at finding the lowest energy cage conformer in a robust way. However, I would expect there are better and more efficient approaches, so editting this protocol is a smart idea.

The `scripts` directory contains many one-off scripts used to produce figures and tables for the manuscript and to setup and analyse DFT validation. But nothing in there is required for the screening.

Contributors and Acknowledgements
---------------------------------

I acknowledge Lukas Turcani for assistance with `stk` development.

License
-------

This project is licensed under the MIT license.
