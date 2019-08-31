#!/usr/bin/python2
# using numpy 1.12
import numpy as np
import sys
import re
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

# minimalDistance(a1, a2, b1, b2)                   -> minimal distance (float)
# octant(p)                                         -> octant           (float)
# isSecondOblateAxis(alpha1, alpha2, beta1, beta2)  -> axis check       (bool)
# ellipsAlignMatrix(a1, a2)                         -> rotation matrix  (3x3 np.matrix)
# isInProlate                                       -> position check   (bool)
# isInOblate                                        -> position check   (bool)

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

def octant(p):
    """ find the octant of the point p"""
    x = p[0]
    y = p[1]
    z = p[2]
    if z > 0:
        if y > 0:
            if x > 0:
                return 1
            else:
                return 2
        else:
            if x > 0:
                return 4
            else:
                return 3
    else:
        if y > 0:
            if x > 0:
                return 5
            else:
                return 6
        else:
            if x > 0:
                return 8
            else:
                return 7

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
    if (abs(DotProdNormed) > maxTors):
   #    print beta1, beta2, "not rectangular, angle = ", np.arccos(DotProdNormed)
        return False
   # print beta1, beta2, "is rectangular."       
    # find nearest point to alpha mid on the potential beta axis by bisection
    # midAlpha = [a2 + 0.5 * dAlph for a2, dAlph in zip(alpha2, dirAlpha)]
    axisDist = minimalDistance(a1, a2, b1, b2)
#    print "Distance of", a1, "<->", a2, "  to  ",  b1, "<->",  b2, "is", axisDist
    #midBeta  = [b2 + 0.5 * dBeta for b2, dBeta in zip(beta2, dirBeta)]
    if axisDist < maxDist:
   #     print b1, "<->",  b2, "is possible axis"
        return True
    else:
   #     print b1, "<->",  b2, "is too far (", axisDist ,") from", a1, "<->", a2, ", maximal allowed distance =", maxDist
        return False

def getEllipsYZRotMatrix(a1, a2):
    """ Creates the rotation matrix to rotate axis a1<->a2 on the x-Axis"""
    adir = a2 - a1
    amid = a1 + 0.5 * adir
    kath = np.sqrt((adir[0] * adir[0]  + adir[1] * adir[1]) / 4.0)
    octantA2 = octant(a2)
    theta = np.arctan( abs( (adir[2]/2) /  kath) )
    #[1, 4, 6, 7 ] => left rotation
    #[2, 3, 5, 8 ] => right rotation
    if octantA2 in [2, 3, 5, 8]: 
        theta = -theta        
    print "theta =" , np.rad2deg(theta)
    RotY = np.matrix( [ [  np.cos(theta), 0.0, np.sin(theta) ],
                        [     0.0       , 1.0,     0.0       ],
                        [ -np.sin(theta), 0.0, np.cos(theta) ]
                      ]) 
    
    psi   = np.arctan( abs( adir[1] / adir[0] ) )
    #[2, 4, 6, 8 ] => left rotation
    #[1, 3, 5, 7 ] => right rotation
    if octantA2 in [1, 3, 5, 7]:
        psi = -psi
    print "psi =" , np.rad2deg(psi)
    RotZ = np.matrix( [ [  np.cos(psi), -np.sin(psi), 0.0 ],
                        [  np.sin(psi),  np.cos(psi), 0.0 ],
                        [      0.0    ,       0.0   , 1.0 ]
                      ])
    return np.asarray( RotY * RotZ )

def getOblateXRotMatrix(aStar1, aStar2):
    """ align the 'pseudo'-orthogonal axis of an oblate
    to the xy-plane by rotation around the x-axis"""
    aStarDir = aStar2 - a1
    aStarmid = aStar1 + 0.5 * aStarDir
    kath = np.sqrt((aStarDir[0] * aStarDir[0]  + aStarDir[1] * aStarDir[1]) / 4.0)
    phi = np.arctan( abs( (aStarDir[2]/2) /  kath) )
    octantAStar2 = octant(aStar2)
    if octantAStar2 in [1, 2, 7, 8]: #
        phi = -phi
    print "phi =" , np.rad2deg(phi)
    RotX = np.matrix( [ [ 1.0,      0.0    ,        0.0 ],
                        [ 0.0,  np.cos(phi), np.sin(phi)],
                        [ 0.0, -np.sin(phi), np.cos(phi)]
                      ])
    return np.asarray( RotX )

#def pointToEllipsAlignment(point, translation, rotation):
#    return rotation * (point - translation)

def isInProlate(sample, alpha, beta):
    """ checks if a sample point (sx,sy,sz) is inside the prolate shape
        with semi=axes alpha > beta. The translation vector and rotation matrix
        have to describe the transformation for aligning the alpha-axis with
        the x of the coordinate system and in setting the center to the origin.
        The fundamental ellipsoidal equation is applied the transformed sample
        point """   
    E    = sample[0] * sample[0] / (alpha * alpha)
    E   += (sample[1] * sample[1] + sample[2] * sample[2] ) / (beta * beta)
    if E > 1.0:
        return False
    else:
        return True
    
def isInOblate(sample, alpha, alphaStar, beta):
#    alphaAv = (alpha + alphaStar) / 2
    alphaAv = alpha
    E = (sample[0] * sample[0]) + (sample[1] * sample[1]) / ( alphaAv * alphaAv )
    E += (sample[2] * sample[2] ) / (beta * beta)
#    print E
    if E > 1.0:
        return False
    else:
        return True

        
#############################
##                         ##
## start main script here  ##
##                         ##
#############################

###################
## settings
###################

np.set_printoptions(precision = 2)
#np.set_printoptions(floatmode = 'fixed') # only available for numpy version>= 1.14

######################################
# parse arguments from command line
######################################

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
elif sys.argv[2] == "--secondAxisTest":
    shapeMode = "secondAxisTest"
else:
    "Please specify a valid shape model (--prolate, --oblate or --all)."

#########################
### Parse *.pdb file
#########################

aAcids = ["ALA", "ARG", "CYS", "GLU", "PHE",
          "GLY", "HIS", "ILE", "LYS", "LEU",
          "MET", "ASN", "PRO", "GLN", "ARG",
          "SER", "THR", "VAL", "TRP", "TYR",
          "PYL", "SEC"]

# reduced set for testing speed-up
#aAcids = ["ALA", "ARG", "CYS", "GLY", "LEU", "Gln", "TYR", "VAL"]

parser = PDBParser(PERMISSIVE=1)
print "Parse ", pdbFile, "..."
#structure = parser.get_structure("1ova", "1ova.pdb")
structure = parser.get_structure(pdbFile[:-4], pdbFile)
model = structure[0]
atomsAll = list(model.get_atoms())

# remove water molecules
atomsAll = list(filter(lambda a: a.get_parent().get_resname() in aAcids, atomsAll))

# Reduce list and keep C-alphas only
atoms1 = list(filter(lambda a: a.get_name() == "CA", atomsAll))
#atoms1 = atomsAll

# print List of search atom
print "\n Start List of atoms:\n"
for i in atoms1:
    print i.get_full_id(), i.get_coord()
print "\n End List of atoms:\n"
    
###############################################
### Find longest connecting axes
### maxDistAtoms1<->maxDistAtoms2
### within the list of atoms
###############################################
atoms2        = atoms1
defaultAtom   = Atom("N", [0.0, 0.0, 0.0], 1.0, 0.0, "d", "d", "d" )
axesNum       = 3
maxDistAtoms1 = [defaultAtom] * axesNum
maxDistAtoms2 = [defaultAtom] * axesNum
maxDists        = [0.0]       * axesNum
minDistOnList   = 0.0
minIndex        = 0
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
          
##########################
# Print search results
##########################
print "\n\n\n Longest distances at: \n"
for i in range(axesNum):
    a1 = maxDistAtoms1[i]
    a2 = maxDistAtoms2[i]
    print "From:" 
#    print  a1.get_parent().get_resname(), " | ", a1.get_name(), " |", a1.get_coord()
    print  a1.get_parent().__repr__(), " | ", a1.get_name(), " |", a1.get_coord()
    print "To:"
    print  a2.get_parent().__repr__(), " | ", a2.get_name(), " |", a2.get_coord()
    print "Distance:"
    print maxDists[i], "\n"
    
#############################
###
###  Prolate Mode
###
#############################

if shapeMode == "Prolate" or shapeMode == "All":
    targetFile = open("testOutProlate.csv","w")        
    print "-----------------------------------"
    print "- Start counting the atoms for     "
    print "- variable beta in prolate         "
    print "-----------------------------------\n"

# Atom counting of prolates outside
    for i in range(axesNum):
        print "+++++++++++++++++++++++++++++++++++++++++"
        print "+++ Check axis a1 <--> a2 \'", i, "\' +++"
        print "+++++++++++++++++++++++++++++++++++++++++\n"
        # take all atoms from the protein
        checkAtoms = [ np.asarray(a.get_coord()) for a in atoms1]
        maxAtoms = len(checkAtoms)
        # take alpha axis      
        a1 = np.asarray(maxDistAtoms1[i].get_coord())
        a2 = np.asarray(maxDistAtoms2[i].get_coord())
        print "a1:", a1, "\na2:", a2, "\namid:", a1 + 0.5 * (a2-a1)
        # shift all atoms from center of a1 <-> a2 axis to origin
        adir = a2 - a1
        alphaCenter = a1 + 0.5 * adir
        shift = -alphaCenter
        atomsShift = [ a + shift for a in checkAtoms ]
        a1Shift = a1 + shift
        a2Shift = a2 + shift
        print "--> Shift amid to origin (and all atoms accordingly) -->"
        print "a1:", a1Shift, "\na2:", a2Shift, "\namid:", a1Shift + 0.5 * (a2Shift-a1Shift)
        # rotate all Shift atoms around origin = new alpha center
        rotMat = getEllipsYZRotMatrix(a1Shift, a2Shift)        
        atomsRot = [ np.matmul(rotMat, a) for a in atomsShift ]
        a1Rot  = np.matmul(rotMat, a1Shift) 
        a2Rot  = np.matmul(rotMat, a2Shift)
        print "@-> Rotate alpha to x-axis @->"
        print "a1:", a1Rot,"\na2:", a2Rot, "\n"
        alpha = abs(a1Rot[0])
        stride = 0.5
        for beta in np.arange(0.5, alpha, stride):
            atomsRot = filter(lambda sample: not isInProlate(sample, alpha, beta), atomsRot)
            if(i == 0 ):
                targetFile.write(str(beta) + "," +  str(len(atomsRot)) + "\n")
    targetFile.close()                            


#############################
###
###  Oblate Mode
###
#############################
if shapeMode == "Oblate" or shapeMode == "All":  
    targetFile = open("testOutOblate.csv","w")        
    print "-----------------------------------"
    print "- Find longest rectangular axes    "
    print "- to first alpha axis              "
    print "-----------------------------------\n"
# Atom counting of prolates outside
    rectAxes1 = [0] * axesNum
    rectAxes2 = [0] * axesNum
    for i in range(axesNum):
        print "+++++++++++++++++++++++++++++++++++++++++"
        print "+++ Check axis a1 <--> a2 \'", i, "\' +++"
        print "+++++++++++++++++++++++++++++++++++++++++\n"
        print "Find rectangular axes by bruteforce. This may take some time..."        
        # take alpha axis      
        a1 = np.asarray(maxDistAtoms1[i].get_coord())
        a2 = np.asarray(maxDistAtoms2[i].get_coord())
        # initialize
        atoms2 = atoms1
        rectAxes1[i]     = [defaultAtom] * axesNum
        rectAxes2[i]     = [defaultAtom] * axesNum
        minIndex         = 0
        minDistOnList    = 0.0
        progressCounter  = 0.0
        maxDists         = [0] * axesNum
        maxAxesDist      = 3.0
        maxAxesTors      = 5.0
        #############################################        
        # go through half matrix of all atom pairs
        # and find longest "almost" rectangular
        # second alpha axes
        #############################################
        for aStar1 in atoms1 :
            atoms2 = atoms2[1:]  #remove a1-equivalent from a2 => check only i > j
            progressCounter += 1
            if (progressCounter % 100 == 0):
                print ("%.2f"% ( progressCounter / len(atoms1) * 100 ) ), "% done..."
            for aStar2 in atoms2:
                # find longest 
                dist = aStar2 - aStar1
                if dist > minDistOnList:
                    if isSecondOblateAxis(a1, a2, aStar1.get_coord(), aStar2.get_coord(), maxAxesDist, maxAxesDist):
                        rectAxes1[i][minIndex] = aStar1
                        rectAxes2[i][minIndex] = aStar2
                        maxDists[minIndex]    = dist
                        minIndex              = maxDists.index(min(maxDists))
                        minDistOnList = maxDists[minIndex]
            
            checkAtoms = [ np.asarray(a.get_coord()) for a in atoms1]
            maxAtoms = len(checkAtoms)            
        for j in range(axesNum):
            print "---------------------------------------------"
            print "--- Check axis a1 <--> a2 \'", i, "\'     ---"
            print "---       and a*1 <--> a*2 \'", j, "\'    ---"
            print "--- Start counting with oblate with       ---"
            print "--- increasing beta axis                  ---"
            print "---------------------------------------------\n"
            aStar1 = np.asarray(rectAxes1[i][j].get_coord())
            aStar2 = np.asarray(rectAxes2[i][j].get_coord())
            # take all atoms from the protein
            checkAtoms = [ np.asarray(a.get_coord()) for a in atoms1]
            maxAtoms = len(checkAtoms)         
            # shift all atoms from center of a1 <-> a2 axis to origin
            adir = a2 - a1
            alphaCenter = a1 + 0.5 * adir
            print "a1:", a1, "\na2:", a2, "\namid:", alphaCenter
            print "a*1:", aStar1,"\na*2:", aStar2, "\n"   
            shift = -alphaCenter
            atomsShift = [ a + shift for a in checkAtoms ]
            a1Shift = a1 + shift
            a2Shift = a2 + shift
            aStar1Shift = aStar1 + shift
            aStar2Shift = aStar2 + shift
            print "--> Shift amid to origin (and all atoms accordingly) -->"            
            print "a1:", a1Shift, "\na2:", a2Shift, "\namid:", a1Shift + 0.5 * (a2Shift-a1Shift)
            print "a*1:", aStar1Shift,"\na*2:", aStar2Shift, "\n"   
            # rotate all Shift atoms around origin = new alpha center
            rotYZMat = getEllipsYZRotMatrix  (a1Shift, a2Shift)        
            atomsRot = [ np.matmul(rotYZMat, a) for a in atomsShift ]
            a1Rot  = np.matmul(rotYZMat, a1Shift) 
            a2Rot  = np.matmul(rotYZMat, a2Shift)
            aStar1Rot  = np.matmul(rotYZMat, aStar1Shift) 
            aStar2Rot  = np.matmul(rotYZMat, aStar2Shift)
            print "@-> Rotate round y- and z-axis to align alpha with x-axis @->"
            print "a1:", a1Rot,"\na2:", a2Rot
            print "a*1:", aStar1Rot,"\na*2:", aStar2Rot, "\n"   
            # rotate all atoms around x to align alpha* with xy-plane
            rotXMat = getOblateXRotMatrix(aStar1, aStar2)
            atomsAligned = [ np.matmul(rotXMat, a) for a in atomsRot ]
            a1Aligned  = np.matmul(rotXMat, a1Rot) 
            a2Aligned  = np.matmul(rotXMat, a2Rot)
            aStar1Aligned  = np.matmul(rotXMat, aStar1Rot) 
            aStar2Aligned  = np.matmul(rotXMat, aStar2Rot)
            print "<-@ Rotate around x-axis to align aStar with <-@"
            print "a1:", a1Aligned,"\na2:", a2Aligned
            print "a*1:", aStar1Aligned,"\na*2:", aStar2Aligned, "\n"   
            alpha = abs( a1Aligned[0] )
            print alpha
            alphaStar = alpha # make length of aDirStar later
#            print a1Aligned, alpha
            stride = 0.5            
            for beta in np.arange(0.5, alpha, stride):
                atomsAligned = filter(lambda sample: not isInOblate(sample, alpha, alphaStar, beta), atomsAligned)
                if (i == 0) and (j == 0):
                    targetFile.write(str(beta) + "," +  str(len(atomsAligned)) + "\n")          
    targetFile.close()
        
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

#######################
# Test cases
#######################
if shapeMode == "secondAxisTest":
    print "-------------------------------------"
    print "- Test rectangular axis recognition: "        
    print "-------------------------------------"
    testA1 = [  0.0,  0.0,  0.0 ]
    testA2 = [ 10.0, 10.0, 10.0 ]
    testB1 = [  0.0,  0.0, 10.1 ]
    testB2 = [ 10.0, 10.0,  0.0 ]
    testMaxDist = 4.0  # angstrom
    testMaxTors = 10.0 # degree
    print "1. test case:", testA1, "<-->", testA2, "    ", testB1, "<-->", testB2
    print  "maximal Distance:", testMaxDist ,"maximal Torsion:", testMaxDist
    print "Distance and torsion small enough?"
    print isSecondOblateAxis(testA1, testA2, testB1, testB2, testMaxDist, testMaxTors), "\n"
    testA1 = [  0.0,  0.0,  0.0 ]
    testA2 = [ 10.0, 10.0, 10.0 ]
    testB1 = [  0.0,  0.0, 10.1 ]
    testB2 = [ 10.0, 10.0,  0.0 ]
    print "2. test case:", testA1, "<-->", testA2, "    ", testB1, "<-->", testB2
    print  "maximal Distance:", testMaxDist ,"maximal Torsion:", testMaxDist
    print "Distance and torsion small enough?"
    print isSecondOblateAxis(testA1, testA2, testB1, testB2, testMaxDist, testMaxTors), "\n"
    testA1 = [  0.0,  0.0,  0.0 ]
    testA2 = [ 10.0, 10.0, 10.0 ]
    testB1 = [  0.0,  0.0, 10.1 ]
    testB2 = [ 10.0, 10.0,  0.0 ]
    print "3. test case:", testA1, "<-->", testA2, "    ", testB1, "<-->", testB2
    print  "maximal Distance:", testMaxDist ,"maximal Torsion:", testMaxDist
    print "Distance and torsion small enough?"
    print isSecondOblateAxis(testA1, testA2, testB1, testB2, testMaxDist, testMaxTors), "\n"
    
