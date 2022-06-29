from pymatgen.core import Structure, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab
from matplotlib import pyplot as plt
import lammps
import lammps

nameOfFile = input("Please type the name of cif file. You can download cif file from "
                   "https://next-gen.materialsproject.org/materials.\n")

baseElement = Structure.from_file(nameOfFile)

print("Please write miller indices of the slab in the form h,k,l.")
hIndex = int(input("h="))

kIndex = int(input("k="))

lIndex = int(input("l="))

mIndices = ([hIndex, kIndex, lIndex])

slabSize = float(input("Please write slab size in Angstrom. Slab size will define total number of atoms in your "
                       "system.\n"))

vacuumSize = float(input("Please write vacuum size in Angstrom.\n"))

slabGen = SlabGenerator(baseElement, mIndices, slabSize, vacuumSize)

createdSlabs = slabGen.get_slabs()

slab = createdSlabs[0]

adsorptionSites = AdsorbateSiteFinder(slab)

adsSites = adsorptionSites.find_adsorption_sites()

noOfAtoms = int(input("Type number of different atoms in adsorbate (# of atoms in reduced formula of adsorbate).\n"))

my_vars = []

for x in range(0, noOfAtoms):
    adsAtoms = input("Type name of atoms.\n")
    my_vars.append(adsAtoms)

print("Write relative positions of atoms in Angstrom.")

axes = ('x=', 'y=', 'z=')


def input_vector(axes):
    values = []
    string = "Enter the coordinate of "
    for coord in axes:
        values.append(float(input(string + coord)))
    return values


vectors = []
for i in range(noOfAtoms):
    vectors.append(input_vector(axes))

adsorbate = Molecule(my_vars, vectors)

dist = float(input("Please write distance from surface of the adsorbate in Angstrom.\n"))

print(adsSites.keys())

chosenAdsSite = input("Please pick adsorption site.\n")

if chosenAdsSite == "ontop":
    x = 1
if chosenAdsSite == "bridge":
    x = 2
if chosenAdsSite == "hollow":
    x = 3
if chosenAdsSite == "all":
    x = 4
adsStructs = adsorptionSites.generate_adsorption_structures(adsorbate, repeat=[1, 1, 1], find_args={"distance": dist})

finalStructure = adsStructs[x - 1]

# Plotting created slabs.

# fig = plt.figure(figsize=[15, 15])

# ax = fig.add_subplot(1, 4, x+1)

# plot = plot_slab(adsStructs[x-1], ax, adsorption_sites=False)
# ax.set_xlim(-10, 10)
# ax.set_ylim(-10, 10)
# plt.show()

h = finalStructure.as_dict()
print(finalStructure)
sites = int(input("Type number of sites written above to inform lammps how many atoms we have in the system.\n"))

coordsString = ""
for x in range(0, sites):
    coords = h["sites"][x]["xyz"]
    atomNames = h["sites"][x]["label"]
    coordsString += atomNames + " 0 "
    for i in range(0, 3):
        coordsString += str(coords[i]) + " "
    coordsString += "\n"

listedCoords = coordsString.splitlines()

nOfTotalAtoms = sites
systemSize = float(input("Please type the system size in Angstrom.\n"))

nOfLammpsTypes = int(input("Please write total number of different atoms in the system.(eg. Ni+CH4 # of different "
                           "atoms = 3)\n")) 
massesOfLammpsSpecies = []
namesOfLammpsSpecies = []

for x in range(0, nOfLammpsTypes):
    lammpsNames = input("Type name of atoms in the lammps simulation.\n")
    namesOfLammpsSpecies.append(lammpsNames)

for x in range(0, nOfLammpsTypes):
    masses = float(input("Type mass of atoms in unit of g/mol\n"))
    massesOfLammpsSpecies.append(masses)

with open("lammpsData.input", "w") as fdata:
    fdata.write('{} and {} input\n'.format(input("Please type the name of base(slab) element for lammps.\n"),
                                           input("Please type the reduced formula of the adsorbate.\n")))
                                                                                                            
    fdata.write('\n')
    fdata.write('         {} atoms\n'.format(sites))
    fdata.write('{} atom types\n'.format(nOfLammpsTypes))
    fdata.write('\n')
    # Specify box dimensions
    fdata.write('{} {} xlo xhi\n'.format(0.0, systemSize))
    fdata.write('{} {} ylo yhi\n'.format(0.0, systemSize))
    fdata.write('{} {} zlo zhi\n'.format(0.0, systemSize))
    fdata.write('\n')
    fdata.write("Masses\n")
    fdata.write('\n')
    for x in range(nOfLammpsTypes):
        fdata.write("  {}     {}   # {}\n".format(x + 1, massesOfLammpsSpecies[x], namesOfLammpsSpecies[x]))
    fdata.write('\n')
    fdata.write('Atoms\n')
    fdata.write('\n')
    for i in range(sites):
        fdata.writelines('{} {}\n'.format(i + 1, listedCoords[i]))

while True:
    cont = input("You need to provide lammps in. and potential data files in source directory. After that open "
                 "lammpsData.input and arrange identifiers and type \"y\".\n")
    if cont == "y":
        break
    else:
        continue

inp = input("Please type the name of the main LAMMPS file.(in.NameOfTheFile)")

lmp = lammps.lammps()

lmp.file(inp)
