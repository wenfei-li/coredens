# Project Name

Coredens

## Description

A simple program which reads valence electron density in Gaussian Cube format, and adds core electron densities from either
- the edflib database (https://github.com/zorkzou/Molden2AIM)
or
- the coredensity library (https://github.com/zezhong-zhang/coredensity)
Then writes the resulting total density to a Gaussian Cube file.

## Table of Contents

- [Usage](#usage)

## Usage

### The main program
Check the directory example/ for how to use this program. The valence density is stored in SPIN1_CHG.cube. I have not made filename a variable, so it has to be named this way.
For nspin = 2 calculations, the up and down densities should be stored as SPIN1_CHG.cube and SPIN2_CHG.cube. The program will detect whether SPIN2_CHG.cube exists in the working directory and set nspin automatically.
In input.dat one must specify which database is used, 0 means edflib and 1 means coredensity library.
If the former is used, no extra input file is needed.
If the latter is used, one must also provide core density files for the elements in the system. Such core density files can be extracted from coredensity library using the Python script script/extract_core_dens.py

### The script extract_core_dens.py
See script/example/input.json for the input file required by the script extract_core_dens.py. Three pieces of information are required: the element, the desired core electron configuration, and the path to the coredensity library.