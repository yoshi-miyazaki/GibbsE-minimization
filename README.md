# GibbsE-minimization
This code is the implementation of Gibbs energy minimization in Miyazaki and Korenaga (ApJ, 2017).

Folder `Gibbs_data`
- Thermodynamic data for species that can be considered in the minimization are stored in this folder.
  - Thermodynamic data for silicates phases are taken from Robie & Hemingway (1995), and those for the rest of the phases are from the JANAF Thermochemical Tables (https://janaf.nist.gov/).

CGgibbsmin.cpp / CGgibbsmin.h
- Class `gibbseminCG` provides the equilibrium composition by Gibbs energy minimization
  - based on Gibbs energy data `MoleculeData_G` and an initial composition
- To calculate the equilibrium composition,
    1. Construct `MoleculeData_G` object from a text file containing list of species
      - `MoleculeData_G   gibbs_list("./spec_atmos.txt", Ts, 1.0);`
    2. Construct `gibbsminCG` object from `MoleculeData_G` and `tensor1d` vector containing molar amount of each species.    
    3. Equilibrium composition is obtained by 
      - `tensor1d<double> nequil = min.getnbest();`

gibbse.cpp / gibbse.h
- Class `Molecule_Gibbs` is an inheritance class of `Molecule`
  - apparent Gibbs energy & phase information at temperature T and pressure P
- Class `MoleculeData_G` contains the list of phase and Gibbs energy


elment.cpp / element.h
- Class defenition of `Element` & `Molecule`, which contains molar mass and other information
  - Make an element or molecule object from `string`

massb.cpp / massb.h
- Class defenition of `massBalance`
  - Construct a mass-balance matrix (see Eq.13 of Miyazaki and Korenaga, 2017)
  - from `initcomp`

initcomp.cpp / initcomp.h
- Class `initcomp` reads the input file including the list of molecules and their molar amount

tensor.h, tensor_s.h, matrix.h, matrix_s.h
- Define `tensor1d`, an array object used in this code

const.h
- List of physical constants
