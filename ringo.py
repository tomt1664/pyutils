#!/usr/bin/python

# ringo.py
# T.Trevethan
# 2015
#
# Script to analyse a covalently bonded structure and output all n-fold
# bonded rings
#
# The structure is supplied as an xyz format file: input.xyz
# 
# The script analyses the input structure and prints the number of 
# 3-fold (triangles) 4-fold (squares) 5-fold (pentagons) 6-fold (hexagons)
# and 7-fold (heptagons)
# 
# The coordinates of the polygons are output to xyz format files if specified 
# at the command line. 
#
# Optional command line arguments:  
#                               
#       -o int n : write n-fold polygon coordinates to nring.xyz
#       -x float : peridic boundary conditions in x-direction
#       -y float : peridic boundary conditions in y-direction
#       -z float : peridic boundary conditions in z-direction
#       -mn float : minimum bond length for connection table (Default 1.0)
#       -mx float : maximum bond length for connection table (Default 1.7) 

import sys
import math 
import argparse
import numpy as np

minb = 1.0
maxb = 1.7

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-o', type = int, help = 'output n-fold ring coordinates')
parser.add_argument('-x', type = float, help = 'x-direction PBC')
parser.add_argument('-y', type = float, help = 'y-direction PBC')
parser.add_argument('-z', type = float, help = 'z-direction PBC')
parser.add_argument('-mn', type = float, help = 'min bond length')
parser.add_argument('-mx', type = float, help = 'max bond length')
variables = parser.parse_args()
nfout = variables.o
xper = variables.x
yper = variables.y
zper = variables.z

if variables.mn:
    minb = variables.mn
if variables.mx:
    maxb = variables.mx

print "Ringo: n-fold ring structure analysis program"

# Open input files
input  = open("input.xyz", 'r')

if nfout:
    if nfout < 3 or nfout > 7:
        print "Error: -o n must be in the range 3 <= n <= 7"
        sys.exit(1)
    print "Set to output "+str(nfout)+" fold rings" 
    outf = str(nfout)+"fold.xyz"
    output = open(outf,'w')

#read input file and store values in array
coords = []

numat = input.readline()
nat = int(numat)
blank = input.readline()

for i in range(nat):
    line = input.readline()
    fdata = line.split()
    data = (float(fdata[1]),float(fdata[2]),float(fdata[3]))
    coords.append(data)

print "Read in "+str(nat)+" atoms"

#create bond table
bond = np.zeros((nat,4),dtype=np.int)

for i in range(nat):
    ib = 0
    for j in range(nat):
        if i != j:
            xdst = coords[i][0] - coords[j][0]
            ydst = coords[i][1] - coords[j][1]
            zdst = coords[i][2] - coords[j][2]
            dst = math.sqrt(xdst*xdst+ydst*ydst+zdst*zdst)
            if dst <= minb:
                print "Error: atoms "+str(i)+" and "+str(j)+" closer than min dist"
                sys.exit(1)
            if dst > minb and dst < maxb:
                if ib > 3:
                    print "Error: atom "+str(i)+" has more than 4 bonds"
                    sys.exit(1)
                bond[i][ib] = j
                ib += 1

#getting the bond list
nn2 = []
num2 = 0
for i in range(nat):
    for ib in range(4):
        if bond[i][ib] > i:
            bdat = (i,bond[i][ib])
            nn2.append(bdat)
            num2 += 1

print "No. 2-body: "+str(num2)

#getting the 3-seq list
nn3 = []
num3 = 0
for i in range(num2):
    for j in range(i+1,num2):
        if nn2[i][0] == nn2[j][0]:
            adat = (nn2[i][1],nn2[i][0],nn2[j][1])
            nn3.append(adat)
            num3 += 1
        elif nn2[i][1] == nn2[j][0]:
            adat = (nn2[i][0],nn2[i][1],nn2[j][1])
            nn3.append(adat)
            num3 += 1
        elif nn2[i][1] == nn2[j][1]:
            adat = (nn2[i][0],nn2[i][1],nn2[j][0])
            nn3.append(adat)
            num3 += 1
        elif nn2[i][0] == nn2[j][1]:
            adat = (nn2[i][1],nn2[i][0],nn2[j][0])
            nn3.append(adat)
            num3 += 1

print "No. 3-body: "+str(num3)

#getting the 4-seq list
nn4 = []
num4 = 0
for i in range(num3):
    for j in range(i+1,num3):
        if nn3[i][1] == nn3[j][0] and nn3[i][2] == nn3[j][1]:
            tdat = (nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][0] == nn3[j][1] and nn3[i][1] == nn3[j][2]:
            tdat = (nn3[j][0],nn3[i][0],nn3[i][1],nn3[i][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][1] == nn3[j][2] and nn3[i][2] == nn3[j][1]:
            tdat = (nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][0])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][0] == nn3[j][1] and nn3[i][1] == nn3[j][0]:
            tdat = (nn3[j][2],nn3[j][1],nn3[j][0],nn3[i][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][1] == nn3[j][0] and nn3[i][0] == nn3[j][1]:
            tdat = (nn3[j][2],nn3[i][0],nn3[i][1],nn3[i][2])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][2] == nn3[j][1] and nn3[i][1] == nn3[j][2]:
            tdat = (nn3[j][0],nn3[i][2],nn3[i][1],nn3[i][0])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][1] == nn3[j][2] and nn3[i][0] == nn3[j][1]:
            tdat = (nn3[j][0],nn3[i][2],nn3[i][1],nn3[i][0])
            nn4.append(tdat)
            num4 += 1
        elif nn3[i][2] == nn3[j][1] and nn3[i][1] == nn3[j][0]:
            tdat = (nn3[i][2],nn3[i][1],nn3[i][0],nn3[j][2])
            nn4.append(tdat)
            num4 += 1

print "No. 4-body: "+str(num4)

#getting the 5-seq list
nn5 = []
num5 = 0
for i in range(num4):
    for j in range(i+1,num4):
        if nn4[i][1] == nn4[j][0] and nn4[i][2] == nn4[j][1] and nn4[i][3] == nn4[j][2]:
            pdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][3])
            nn5.append(pdat)
            num5 += 1
        elif nn4[i][0] == nn4[j][1] and nn4[i][1] == nn4[j][2] and nn4[i][2] == nn4[j][3]:
            pdat = (nn4[j][0],nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3])
            nn5.append(pdat)
            num5 += 1
        elif nn4[i][1] == nn4[j][3] and nn4[i][2] == nn4[j][2] and nn4[i][3] == nn4[j][1]:
            pdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][0])
            nn5.append(pdat)
            num5 += 1
        elif nn4[i][0] == nn4[j][2] and nn4[i][1] == nn4[j][1] and nn4[i][2] == nn4[j][0]:
            pdat = (nn4[j][3],nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3])
            nn5.append(pdat)
            num5 += 1

print "No. 5-body: "+str(num5)

#find all the triangles
r3 = []
numr3 = 0
for i in range(num4):
    if nn4[i][0] == nn4[i][3]:
        tdat = (nn4[i][0],nn4[i][1],nn4[i][2])
        r3.append(tdat)
        numr3 += 1

print "No. 3-fold rings: "+str(numr3/3)

#find all the squares
r4 = []
numr4 = 0
for i in range(num3):
    for j in range(i+1,num3):
        if nn3[i][0] == nn3[j][0] and nn3[i][2] == nn3[j][2]:
            sdat = (nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][1])
            r4.append(sdat)
            numr4 += 1
        elif nn3[i][0] == nn3[j][2] and nn3[i][2] == nn3[j][0]:
            sdat = (nn3[i][0],nn3[i][1],nn3[i][2],nn3[j][1])
            r4.append(sdat)
            numr4 += 1

print "No. 4-fold rings: "+str(numr4/4)

#find all the pentagons
r5 = []
numr5 = 0
for i in range(num4):
    for j in range(num3):
        if nn4[i][0] == nn3[j][0] and nn4[i][3] == nn3[j][2]:
            if nn4[i][1] != nn3[j][1] and nn4[i][2] != nn3[j][1]:
                sdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn3[j][1])
                r5.append(sdat)
                numr5 += 1
        elif nn4[i][3] == nn3[j][0] and nn4[i][0] == nn3[j][2]: 
            if nn4[i][1] != nn3[j][1] and nn4[i][2] != nn3[j][1]:
                sdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn3[j][1])
                r5.append(sdat)
                numr5 += 1

print "No. 5-fold rings: "+str(numr5/5)

#find all the hexagons
r6 = []
numr6 = 0
for i in range(num4):
    for j in range(i+1,num4):
        if nn4[i][0] == nn4[j][0] and nn4[i][3] == nn4[j][3]:
            if nn4[i][1] != nn4[j][1] and nn4[i][2] != nn4[j][2]:
                hdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][2],nn4[j][1])
                r6.append(hdat)
                numr6 += 1
        elif nn4[i][0] == nn4[j][3] and nn4[i][3] == nn4[j][0]:
            if nn4[i][2] != nn4[j][1] and nn4[i][1] != nn4[j][2]:
                hdat = (nn4[i][0],nn4[i][1],nn4[i][2],nn4[i][3],nn4[j][1],nn4[j][2])
                r6.append(hdat)
                numr6 += 1

print "No. 6-fold rings: "+str(numr6/6)

#find all the heptagons
r7 = []
numr7 = 0
for i in range(num5):
    for j in range(num4):
        if nn5[i][0] == nn4[j][0] and nn5[i][4] == nn4[j][3]:
            if nn5[i][1] != nn4[j][1] and nn5[i][3] != nn4[j][2]:
                hdat = (nn5[i][0],nn5[i][1],nn5[i][2],nn5[i][3],nn5[i][4],nn4[j][2],nn4[j][1])
                r7.append(hdat)
                numr7 += 1
        elif nn5[i][0] == nn4[j][3] and nn5[i][4] == nn4[j][0]:
            if nn5[i][1] != nn4[j][2] and nn5[i][3] != nn4[j][1]:
                hdat = (nn5[i][0],nn5[i][1],nn5[i][2],nn5[i][3],nn5[i][4],nn4[j][1],nn4[j][2])
                r7.append(hdat)
                numr7 += 1

print "No. 7-fold rings: "+str(numr7/7)

#write selected ring structures to output file
if nfout == 3:
    outline = str(numr3*3)+"\n"
    output.write(outline)
    output.write(" \n")
    for i in range(numr3):
        for j in range(3):
            outline = "3 "+str(coords[r3[i][j]][0])+" "+str(coords[r3[i][j]][1])+" "+str(coords[r3[i][j]][2])+"\n"
            output.write(outline)
    print "Output: 3-fold rings"
if nfout == 4:
    outline = str(numr4*4)+"\n"
    output.write(outline)
    output.write(" \n")
    for i in range(numr4):
        for j in range(4):
            outline = "4 "+str(coords[r4[i][j]][0])+" "+str(coords[r4[i][j]][1])+" "+str(coords[r4[i][j]][2])+"\n"
            output.write(outline)
    print "Output: 4-fold rings"
elif nfout == 5:
    outline = str(numr5*5)+"\n"
    output.write(outline)
    output.write(" \n")
    for i in range(numr5):
        for j in range(5):
            outline = "5 "+str(coords[r5[i][j]][0])+" "+str(coords[r5[i][j]][1])+" "+str(coords[r5[i][j]][2])+"\n"
            output.write(outline)
    print "Output: 5-fold rings"
elif nfout == 6:
    outline = str(numr6*6)+"\n"
    output.write(outline)
    output.write(" \n")
    for i in range(numr6):
        for j in range(6):
            outline = "6 "+str(coords[r6[i][j]][0])+" "+str(coords[r6[i][j]][1])+" "+str(coords[r6[i][j]][2])+"\n"
            output.write(outline)
    print "Output: 6-fold rings"
elif nfout == 7:
    outline = str(numr7*7)+"\n"
    output.write(outline)
    output.write(" \n")
    for i in range(numr7):
        for j in range(7):
            outline = "6 "+str(coords[r7[i][j]][0])+" "+str(coords[r7[i][j]][1])+" "+str(coords[r7[i][j]][2])+"\n"
            output.write(outline)
    print "Output: 7-fold rings"
