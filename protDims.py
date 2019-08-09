#!/usr/bin/python2
import numpy as np
import sys
import re
#from scipy.optimize import minimize
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

def distance(point1, point2):
    vx = point2[0] - point1[0]
    vy = point2[1] - point1[1]
    vz = point2[2] - point1[2]
    return sqrt( vx**2 + vy**2 + vz**2 )

def minimalDistance(a1, a2, b1, b2):
    """Get the distance between nearest points to each other 
     on to axes a and b without necesary intersection 
     under the assumption that the points are located\
     between a1,a2 ; b1,b respectively """
    adir = a2 - a1
    bdir = b2 - b1
    amid = a1 + 0.5 * adir
    s    = b1 - amid
    A    = np.dot(bdir, bdir)
    B_2  = np.dot(bdir,    s)
    lambda_beta = - B_2 / A
    bOpt = lambda_beta * bdir + b1
    s = a1 - bOpt
    A    = np.dot(adir, adir)
    B_2  = np.dot(adir, s)
    lambda_alpha = - B_2 / A
    aOpt = lambda_alpha * adir + a1
    Delta = bOpt - aOpt
    return np.sqrt(np.dot(Delta, Delta))

def isProlateBetaAxis(alpha1, alpha2, beta1, beta2, maxDist, maxTorsAngle):
    """Check if an axis from beta1 to beta2 is nearly perpendicular with a maximal 
    distance to an axis from alpha1 to alpha2 and if their distance is under maxDist"""       
    a1 = np.asarray(alpha1)
    a2 = np.asarray(alpha2)
    b1 = np.asarray(beta1)
    b2 = np.asarray(beta2)
    #lent = alpha1 - beta1
    #alphaDotBet = dirAlpha[0] * dirBeta[0] + dirAlpha[1] * dirBeta[1] + dirAlpha[2] * dirBeta[2]
    alphaDotBet = 0.01
    maxTors = np.cos(maxTorsAngle * 3.14 / 180.0)
    if abs(alphaDotBet) > maxTors:
       print "Not rectangular."
       return False

    print "Rectangular; Test maximal distance :"       
   # find nearest point to alpha mid on the potential beta axis by bisection
   # midAlpha = [a2 + 0.5 * dAlph for a2, dAlph in zip(alpha2, dirAlpha)]
    minDist = minimalDistance(a1, a2, b1, b2)
    print minDist
    #midBeta  = [b2 + 0.5 * dBeta for b2, dBeta in zip(beta2, dirBeta)]


#############################
##                         ##
## start main script here  ##
##                         ##
#############################

# parse arguments
if len(sys.argv) != 3:# or sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print "Usage of protDims:"
    print ">>> $~ ./protDims.py <myStructure.pdb> <shapeModel>"
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

# reduced for testing:

aAcids = ["ALA", "ARG", "CYS", "GLY", "LEU", "Gln", "TYR", "VAL"]

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
#C alpha test:
atoms1 = list(filter(lambda a: a.get_name() == "CA", atoms1))
print len(atoms1)
atoms2 = atoms1
#from Bio.PDB import Atom
defaultAtom = Atom("N", [0.0, 0.0, 0.0], 1.0, 0.0, "d", "d", "d" )
axesNum = 5
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
    a2 = maxDistAtoms2[i]
    print "From:" 
    print a1.get_parent().get_resname(), " | ", a1.get_name(), " |", a1.get_coord()
    print "To:"
    print a2.get_parent().get_resname(), " | ", a2.get_name(), " |", a2.get_coord()
    print "Distance:"
    print maxDists[i]
# get axis centers

# centers for each axis i => i X [x,y,z] 
centers = [[ 0.5 * (x1+x2) for x1, x2 in zip(a1.get_coord(), a2.get_coord()) ]
           for a1, a2 in zip(maxDistAtoms1, maxDistAtoms2) ]

# test call
print "TestProlate", isProlateBetaAxis( [0, 0, 0], [10,10,10], [1,1,1], [11, 11, 11], 4.0, 10.0 )


if shapeMode == "Prolate" or shapeMode == "All":
#    maxDistAtoms1[0]
#    maxDistAtoms2[0]
    atoms2 = atoms1

    print "Find orthogonal axes"   
    for a1 in atoms1:
        if (atomCounter % 100 == 0) :
            print ("%.2f"% ( atomCounter / len(atoms1) * 100 ) ), "% done..."
    atoms2 = atoms2[1:] #remove a1-equivalent from a2 => check only i > j    
    #     
if shapeMode == "Oblate"  or shapeMode == "All": pass
    # define "search cylinder" around beta-axis
#define 
