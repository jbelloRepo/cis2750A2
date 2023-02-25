# import molecule
from molecule import molecule, atom, bond
import math

# Define constants
radius = {'H': 25, 'C': 40, 'O': 40, 'N': 40}
element_name = {'H': 'grey', 'C': 'black', 'O': 'red', 'N': 'blue'}
header = """<svg version="1.1" width="1000" height="1000"xmlns="http://www.w3.org/2000/svg">"""
footer = """</svg>"""
offsetx = 500
offsety = 500


class Atom:  # Create an ATOM class
    def __init__(self, atom):
        self.atom = atom
        self.z = atom.z

    def __str__(self):  # Define the string
        return f"element: {self.atom.element} x: {self.atom.x} y: {self.atom.y} z:{self.atom.z}"

    def svg(self):  # svg method: description of circle
        cx = self.atom.x * 100.0 + offsetx
        cy = self.atom.y * 100.0 + offsety
        r = radius[self.atom.element]
        fill = element_name[self.atom.element]
        return f"  <circle cx=\"{cx:.2f}\" cy=\"{cy:.2f}\" r=\"{r}\" fill=\"{fill}\"/>\n"


class Bond:
    def __init__(self, bond):
        self.bond = bond
        self.z = bond.z

    def __str__(self):  # Define the string
        return f"a1: {self.bond.a1} a2: {self.bond.a2} x1: {self.bond.x1} x2: {self.bond.x2} y1: {self.bond.y1} y2: {self.bond.y2} z:{self.bond.z}"

    def svg(self):
        x1 = self.bond.x1 * 100.0 + offsetx
        y1 = self.bond.y1 * 100.0 + offsety
        x2 = self.bond.x2 * 100.0 + offsetx
        y2 = self.bond.y2 * 100.0 + offsety
        dx = x2 - x1
        dy = y2 - y1
        length = math.sqrt(dx*dx + dy*dy)
        ux = dx / length
        uy = dy / length
        x3 = x1 + 10*uy
        y3 = y1 - 10*ux
        x4 = x1 - 10*uy
        y4 = y1 + 10*ux
        x5 = x2 + 10*uy
        y5 = y2 - 10*ux
        x6 = x2 - 10*uy
        y6 = y2 + 10*ux
        return f'  <polygon points="{x3:.2f},{y3:.2f} {x4:.2f},{y4:.2f} {x6:.2f},{y6:.2f} {x5:.2f},{y5:.2f}" fill="green"/>\n'


class Molecule(molecule):
    def __init__(self, molecule):
        self.molecule = molecule
        super().__init__()
        # self.name = name

    def __len__(self):
        return len(self.list)

    # def get_name(self):
    #     return self.name

    def __str__(self):
        return f"atom max: {self.atom_max} atom no: {self.atom_no} bond max: {self.bond_max} bond no: {self.bond_no}"

    # def __str__(self):  # Define the string for debugging
    #     res = "Atoms:\n"
    #     for atom in self.atoms:
    #         res += str(atom) + "\n"
    #     res += "\nBonds:\n"
    #     for bond in self.bonds:
    #         res += str(bond) + "\n"
    #     return res

    def svg(self):
        atoms = self.atoms.copy()
        bonds = self.bonds.copy()
        atoms.sort(key=lambda x: x.z)
        bonds.sort(key=lambda x: min(x.a1.z, x.a2.z))
        svg_str = header
        while atoms or bonds:
            if not atoms:
                svg_str += bonds.pop(0).svg()
            elif not bonds:
                svg_str += atoms.pop(0).svg()
            elif atoms[0].z < min(bonds[0].a1.z, bonds[0].a2.z):
                svg_str += atoms.pop(0).svg()
            else:
                svg_str += bonds.pop(0).svg()
            svg_str += footer
            return svg_str

    @classmethod
    def parse(self, fileobj):
        self = molecule()
        # self.append_atom("O", 2.5369, -0.1550, 0.0000)
        # atom = self.get_atom(1)
        # newAtom = Atom(atom)
        print(self.atom_no)
        atomCount, bondCount = None, None  # declare as object to make Iterable
        lineCount = 0
        for line in fileobj:
            lineCount += 1
            if lineCount == 4:  # read after the header file
                atomCount = int(line[:3])
                bondCount = int(line[3:6])
                # print(f"atom count:{atomCount}  bond count:{bondCount}")
            elif atomCount is not None and self.atom_no < atomCount and (lineCount < (4+atomCount+1)):
                x, y, z, element = line.split()[:4]
                x = float(x)
                y = float(y)
                z = float(z)
                self.append_atom(element, x, y, z)
            elif bondCount is not None and self.bond_no < bondCount and (lineCount > (4+atomCount)) and (lineCount < (4+atomCount+1+bondCount)):
                bondVar = (int(line[:3].strip()), int(line[3:6].strip()), int(line[6:9].strip()))
                self.append_bond(bondVar[0], bondVar[1], bondVar[2])
                # print(bondVar)
        
        print(self.atom_no)
        print(self.bond_no)
                
                # break

        # for line in fileobj:
        #     print(line.rstrip('\n'))

        # mol = Molecule(molecule)
        # atom_count = int(fileobj.readline().strip())
        # bond_count = int(fileobj.readline().strip())
        # fileobj.readline()  # ignore next line

        # # parse atoms
        # for i in range(atom_count):
        #     line = fileobj.readline()
        #     x, y, z, sym = line[:10], line[10:20], line[20:30], line[31:33]
        #     atom = Atom(float(x), float(y), float(z), sym)
        #     mol.append_atom(atom)

        # # parse bonds
        # for i in range(bond_count):
        #     line = fileobj.readline()
        #     a1_idx, a2_idx, bond_type = line[:3], line[3:6], line[6:9]
        #     a1 = mol.atoms[int(a1_idx) - 1]
        #     a2 = mol.atoms[int(a2_idx) - 1]
        #     bond = Bond(a1, a2, int(bond_type))
        #     mol.append_bond(bond)

        # return mol


# mol = molecule()  # create a new molecule object
# mol.append_atom("O", 2.5369, -0.1550, 0.0000)
# mol.append_atom("H", 3.0739, 0.1550, 0.0000)
# mol.append_atom("H", 2.0000, 0.1550, 0.0000)

# mol.append_bond(1, 2, 1)
# mol.append_bond(1, 3, 1)

# atom = mol.get_atom(1)
# bond = mol.get_bond(1)
# z = atom.z


# newAtom = Atom(atom)
# newBond = Bond(bond)
# print(newAtom.svg())
# print(newBond.svg())
# # print()

# from molecule import Molecule, Atom, Bond

# Create a new molecule
mol1 = Molecule(molecule)
# mol1.append_atom("O", 2.5369, -0.1550, 0.0000)
# mol1.append_atom("H", 3.0739, 0.1550, 0.0000)
# mol1.append_atom("H", 2.0000, 0.1550, 0.0000)

# mol1.append_bond(1, 2, 1)
# mol1.append_bond(1, 3, 1)

with open('caffeine-3D-structure-CT1001987571.sdf', 'r') as sdfile:
    mol1.parse(sdfile)
    print(mol1)

# print molecule information

