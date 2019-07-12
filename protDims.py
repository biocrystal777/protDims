#
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

aAcids = ["ALA", "CYS", "GLU", "PHE", "GLY", "HIS", "ILE",
          "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
          "SER", "THR", "VAL", "TRP", "TYR" ]

parser = PDBParser(PERMISSIVE=1)
#builder = StructureBuilder()
structure = parser.get_structure("1ova", "1ova.pdb")
model  = structure[0]
atoms1 = list(model.get_atoms())
## filter water molecules
print len(atoms1)
atoms1  = list(filter(lambda a: a.get_parent().get_resname() in aAcids, atoms1))
print len(atoms1)
atoms2 = atoms1
#from Bio.PDB import Atom
defaultAtom = Atom("N", [0.0, 0.0, 0.0], 1.0, 0.0, "d", "d", "d" )
axesNum = 20
maxDistAtoms1 = [defaultAtom] * axesNum
maxDistAtoms2 = [defaultAtom] * axesNum
maxDists =      [0.0]         * axesNum
# for i in range(axesNum) :
# maxDistAtomPairs.append(defaultAtom, defaultAtom, 0.0)

#maxDists = [0.0] * axesNum

# find longest alpha axes

minDistOnList = 0.0
minIndex      = 0
for a1 in atoms1 :
    atoms2 = atoms2[1:]  #remove a1-equivalent from a2 => check only i > j    
#    print "Test distances with |", a1.get_parent().get_resname(),  "|", a1.get_name(), "|", a1.get_coord()
    for a2 in atoms2 :
       dist = a2 - a1
       if dist > minDistOnList :
          maxDistAtoms1[minIndex] = a1
          maxDistAtoms2[minIndex] = a2
          maxDists[minIndex]      = dist
          minIndex                = maxDists.index(min(maxDists))
          minDistOnList = maxDists[minIndex]
print "\n\n\n Longest distances at:"
for i in range(axesNum):
    a1 = maxDistAtoms1[i]
    a2 = maxDistAtoms1[i]
    print "From:" 
    print a1.get_parent().get_resname(), " | ", a1.get_name(), " |", a1.get_coord()
    print "To:"
    print a2.get_parent().get_resname(), " | ", a2.get_name(), " |", a2.get_coord()
    print "Distance:"
    print maxDists[i]
    print distance
    # get axis midpoints


    
