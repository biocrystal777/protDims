#!/usr/bin/python2
import numpy as np
import sys
import re
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

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

def isSecondOblateAxis(alpha1, alpha2, beta1, beta2, maxDist, maxTorsAngle):
    """Check if an axis from beta1 to beta2 is nearly perpendicular with a maximal 
    distance to an axis from alpha1 to alpha2 and if their distance is under maxDist"""       
    a1 = np.asarray(alpha1)
    a2 = np.asarray(alpha2)
    b1 = np.asarray(beta1)
    b2 = np.asarray(beta2)
    #lent = alpha1 - beta1
    adir = a2 - a1
    bdir = b2 - b1
    aLength = np.sqrt ( np.dot(adir, adir) )
    bLength = np.sqrt ( np.dot(bdir, bdir) )
    DotProdNormed = np.dot(adir, bdir) / ( aLength * bLength )   
    maxTors = np.cos( np.radians( maxTorsAngle ))
    if abs(DotProdNormed) > maxTors:
       print beta1, beta2, "not rectangular, angle = ", np.arccos(DotProdNormed)
       return False
    print beta1, beta2, "is rectangular."       
    # find nearest point to alpha mid on the potential beta axis by bisection
    # midAlpha = [a2 + 0.5 * dAlph for a2, dAlph in zip(alpha2, dirAlpha)]
    axisDist = minimalDistance(a1, a2, b1, b2)
    print "Distance of", a1, "<->", a2, "  to  ",  b1, "<->",  b2, "is", axisDist
    #midBeta  = [b2 + 0.5 * dBeta for b2, dBeta in zip(beta2, dirBeta)]
    if (axisDist < maxDist):
        print b1, "<->",  b2, "is possible axis"
        return True
    else:
        print b1, "<->",  b2, "is too far (", axisDist ,") from", a1, "<->", a2, ", maximal allowed distance =", maxDist
        return False

def ellipsAlignMatrix(a1, a2):
    adir = a2 - a1
    amid = a1 + 0.5 * adir
    kath = np.sqrt((adir[0] * adir[0]  + adir[1] * adir[1]) / 4.0)
    theta = -np.arctan( k1 / ( adir[2]/2)  )
    RotY = np.matrix( [ [  np.cos(theta), 0.0, np.sin(theta) ],
                        [     0.0       , 1.0,     0.0       ],
                        [ -np.sin(theta), 0.0, np.cos(theta) ]
                      ]) 
    psi   = -np.arctan( adir[1] / adir[0]  )
    RotZ = np.matrix( [ [  np.cos(psi), -np.sin(psi), 0.0 ],
                        [  np.sin(psi),  np.cos(psi), 0.0 ],
                        [      0.0    ,       0.0   , 1.0 ]
                      ])  
    return RotY * RotZ

def pointToEllipsAlignment(point, translation, rotation):
    return rotation * (point - translation)

def isWithinProlate(a1, a2, betaLength, normedSampl): pass
#########################################################
# transform coordinate system for sample in such way
# that the center a1<->a2 is on (0,0,0)
# and the a1<->a2 is aligned with x axis
#########################################################

    

#############################
##                         ##
## start main script here  ##
##                         ##
#############################

#####################
# parse arguments from command line
#####################

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

####################
### Parse *.pdb file
####################

aAcids = ["ALA", "ARG", "CYS", "GLU", "PHE",
          "GLY", "HIS", "ILE", "LYS", "LEU",
          "MET", "ASN", "PRO", "GLN", "ARG",
          "SER", "THR", "VAL", "TRP", "TYR",
          "PYL", "SEC"]

# reduced set for testing speed-up
#aAcids = ["ALA", "ARG", "CYS", "GLY", "LEU", "Gln", "TYR", "VAL"]

parser = PDBParser(PERMISSIVE=1)
print "Parse ", pdbFile, "..."
structure = parser.get_structure("1ova", "1ova.pdb")
model = structure[0]
atomsAll = list(model.get_atoms())

# remove water molecules
atomsAll  = list(filter(lambda a: a.get_parent().get_resname() in aAcids, atomsAll))

# Reduce list and keep C-alphas only
atoms1 = list(filter(lambda a: a.get_name() == "CA", atomsAll))

###############################################
### Find longest connecting axes
### maxDistAtoms1<->maxDistAtoms2
### within the list of atoms
###############################################
atoms2 = atoms1
defaultAtom = Atom("N", [0.0, 0.0, 0.0], 1.0, 0.0, "d", "d", "d" )
axesNum = 5
maxDistAtoms1 = [defaultAtom] * axesNum
maxDistAtoms2 = [defaultAtom] * axesNum
maxDists      = [0.0]         * axesNum
minDistOnList = 0.0
minIndex      = 0
progressCounter = 0.0
print "Find longest axes by raw brute force. This will take some time."
for a1 in atoms1 :
    atoms2 = atoms2[1:]  #remove a1-equivalent from a2 => check only i > j    
    #    print "Test distances with |", a1.get_parent().get_resname(),  "|", a1.get_name(), "|", a1.get_coord()
    progressCounter += 1
    if (progressCounter % 100 == 0):
        print ("%.2f"% ( progressCounter / len(atoms1) * 100 ) ), "% done..."
    for a2 in atoms2 :
       dist = a2 - a1
       if dist > minDistOnList :
          maxDistAtoms1[minIndex] = a1
          maxDistAtoms2[minIndex] = a2
          maxDists[minIndex]      = dist
          minIndex                = maxDists.index(min(maxDists))
          minDistOnList = maxDists[minIndex]
print "\n\n\n Longest distances at: \n"
for i in range(axesNum):
    a1 = maxDistAtoms1[i]
    a2 = maxDistAtoms2[i]
    print "From:" 
    print a1.get_parent().get_resname(), " | ", a1.get_name(), " |", a1.get_coord()
    print "To:"
    print a2.get_parent().get_resname(), " | ", a2.get_name(), " |", a2.get_coord()
    print "Distance:"
    print maxDists[i]

if shapeMode == "Prolate" or shapeMode == "All":
#############################
# Atom counting on 
#############################
    for i in range(axesNum):
        adir = maxDistAtoms2[i] - maxDistAtoms1[i]
        
        
#    maxDistAtoms1[0]
#    maxDistAtoms2[0]
#    atoms2 = atoms1
#
#    print "Find orthogonal axes"   
#    for a1 in atoms1:
#        if (progressCounter % 100 == 0) :
#            print ("%.2f"% ( progressCounter / len(atoms1) * 100 ) ), "% done..."
#    atoms2 = atoms2[1:] #remove a1-equivalent from a2 => check only i > j    
    #     
if shapeMode == "Oblate"  or shapeMode == "All":
    
# test call for recognition of orthogonal beta axis
    print "-------------------------------------"
    print "\nTest Axis recognition:", isSecondOblateAxis( [0.0, 0.0, 0.0], [10.0,10.0,10.0], [0.0, 0.0, 10.1], [10, 10, 0], 4.0, 10.0 ), "\n"
    print "-------------------------------------"
