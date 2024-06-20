# Based on the core electron configuration of an element, accumulate
# the core electron density from corresponding core orbitals.

import numpy as np
import h5py as h5
import json as js
from scipy import integrate, interpolate
import matplotlib.pyplot as plt

ELEMENT_NAMES = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]

# Read the core electron configuration from a JSON file
with open('input.json') as input_file:
    input_data = js.load(input_file)

element = input_data['element']
zat = ELEMENT_NAMES.index(element)+1
core_config = input_data['core_config']
database_dir = input_data['database_path']

core_orb_list = [ core_config[i:i+2] for i in range(0, len(core_config), 2) ]

# Read orbital densities from the database and accumulate to get tot core density
with h5.File(database_dir + '/' + element + '.h5', 'r') as f:
    indices = list(f.keys())
    indices.sort(key=lambda x: int(x))

    r0 = np.asarray(f[indices[-1]]['r'][:])
    core_dens = np.zeros_like(r0)

    for index in indices:
        orbital = f[index].attrs['orbital']
        if(any(orb in orbital for orb in core_orb_list)):
            r = np.asarray(f[index]['r'][:])
            density = f[index]['density'][:]
            tck = interpolate.splrep(r, density)
            density_spline = interpolate.splev(r0, tck)
            core_dens = core_dens + density_spline

    nelec = integrate.simps(core_dens, r0)
    print("Number of core electrons is : "+str(nelec)+", make sure this is as expected.")

    filename = 'core_dens_'+str(zat)+'.dat'
    np.savetxt(filename, np.column_stack((r0, core_dens)))