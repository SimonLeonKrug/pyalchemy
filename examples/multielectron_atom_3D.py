import numpy as np
from pyalchemy import kernel_3D, partial_v_mol_3D
from numba import prange # To speed up the summation
from pyscf import gto, scf, qmmm, dft # To get the electron density and energies
import basis_set_exchange as bse # For easy access to different basis sets


# ---------------------------------Parameters-----------------------------------

# number of electrons
n_e = 8 #must be an integer!

# nuclear charges
Z_A = 8
Z_B = 7

basis_set = 'def2-TZVP'
reference_atom = 'Xe'


# ------------------------------External potential------------------------------

# External potential of the initial system
def partial_v_A(n_x, n_y, n_z, x, y, z):
    return partial_v_mol_3D([[Z_A,0,0,0]], n_x, n_y, n_z, x, y, z)

# External potential of the final system
def partial_v_B(n_x, n_y, n_z, x, y, z):
    return partial_v_mol_3D([[Z_B,0,0,0]], n_x, n_y, n_z, x, y, z)


# --------------------------Multi-electron calculation--------------------------

initial_charge = int(Z_A - n_e)

elements = {'H':1, 'He':2,
'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10,
'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18,
'K':19, 'Ca':20, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30,
'Rb':37, 'Sr':38, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53,'Xe':54,
'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44,'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
'Cs':55, 'Ba':56, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86,
'La':57,
'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80}

inv_elements = {v: k for k, v in elements.items()}

# Get the necessary basis functions
basis_specs = {}
for i in range(1,86+1):
    try:
        basis_specs.update({inv_elements[i]: gto.basis.load(bse.get_basis(basis_set, fmt='nwchem'), reference_atom)})
    except:
        pass

#Monkey patching PySCF's qmmm
def add_qmmm(calc, molecule, Z):
    mf_inter = qmmm.mm_charge(calc, molecule.atom_coords(), Z)
    def energy_nuc(self):
        q = molecule.atom_charges().copy().astype(np.float64)
        q += Z
        return molecule.energy_nuc(q)
    mf_inter.energy_nuc = energy_nuc.__get__(mf_inter, mf_inter.__class__)
    return mf_inter

# Initialize the initial atom
mol = gto.M(atom = str(int(Z_A))+' 0 0 0',
            unit='Bohr',
            charge=initial_charge,
            spin=n_e%2,
            symmetry = True,
            verbose=0,
            basis = basis_specs)

# Start an unrestricted Hartree-Fock calculation and stack the non-integer charges
# on top of the nucleus
mf = add_qmmm(scf.UHF(mol), mol, [Z_A - int(Z_A)])
mf.kernel()

# Get the density matrix
dm = sum(mf.make_rdm1())

# Wrap initial density in a callable
def rho_initial(x,y,z):
    ao_value = dft.numint.eval_ao(mol, [[x,y,z]], deriv=0)
    return dft.numint.eval_rho(mol, ao_value, dm, xctype='LDA')[0]

# Find the energy difference of the SCF calculation for comparison
def Delta_E_SCF():
    energy_initial = mf.e_tot
    # Stack the difference in charges into the nucleus
    mf_final = add_qmmm(scf.UHF(mol), mol, [Z_B - int(Z_A)])
    mf_final.kernel()
    energy_final = mf_final.e_tot
    return energy_final - energy_initial


# ----------------------------Integration grid----------------------------------

#Generate grid points and weights with Becke-Lebedev scheme
grid = dft.gen_grid.Grids(mol)
grid.level = 3 # 3 suffices, 5 is precise
grid.build()


# -------------------------Alchemical Integral Transform------------------------

# the true value from the SCF solution
print(Delta_E_SCF())

# the predicted value from AIT
sum = 0
for c in prange(len(grid.coords)):
    sum += grid.weights[c]*rho_initial(*grid.coords[c])*kernel_3D(partial_v_A, partial_v_B, *grid.coords[c])
print(sum)
