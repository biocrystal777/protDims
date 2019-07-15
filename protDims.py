#!/usr/bin/python2
import sys
import re
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

# parse arguments
if len(sys.argv) != 3:# or sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print "Usage of protDims:"
    print ">>> $~ ./protDims.py  <myStructure.pdb> <shapeModel>"
    print "A <shapeModel> has to be provided as either --prolate, --oblate or --all "
    print "Example: "
    print ">>> $~ ./protDims  1ova.pdb"
    print "(Unfortunately) this script is written with Python 2.7"
    sys.exit()
    
print sys.argv[1]
if re.match( "[a-zA-Z0-9_]*.pdb", sys.argv[1]):
    pdbFile = sys.argv[1]
else:
    print sys.argv[1], "does not seem to be a *.pdb file."
    sys.exit()
    
if  sys.argv[2] == "--oblate":
    shapeMode = "Oblate"
    print "Run with oblate shape model..."
elif sys.argv[2] == "--prolate":
    shapeMode = "Prolate"
    print "Run with prolate shape model..."
elif sys.argv[2] == "--all":
    shapeMode =  "All"
    print "Run with both shape models..."
else:
    "Please specify a valid shape model (--prolate, --oblate or --all)."

aAcids = ["ALA", "ARG", "CYS", "GLU", "PHE",
          "GLY", "HIS", "ILE", "LYS", "LEU",
          "MET", "ASN", "PRO", "GLN", "ARG",
          "SER", "THR", "VAL", "TRP", "TYR",
          "PYL", "SEC"]

parser = PDBParser(PERMISSIVE=1)
#builder = StructureBuilder()
print "Parse ", pdbFile, "..."
structure = parser.get_structure("1ova", "1ova.pdb")
model = structure[0]
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
atomCounter   = 0.0 
print "Find longest axes by raw brute force. This will take some time."
for a1 in atoms1 :
    atoms2 = atoms2[1:]  #remove a1-equivalent from a2 => check only i > j    
    #    print "Test distances with |", a1.get_parent().get_resname(),  "|", a1.get_name(), "|", a1.get_coord()
    atomCounter += 1.0
    if (atomCounter % 100 == 0) :
        print ("%.2f"% ( atomCounter / len(atoms1) * 100 ) ), "% done..."
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
# get axis midpoints

    
