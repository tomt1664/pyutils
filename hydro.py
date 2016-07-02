#!/usr/bin/python

# hydro.py
# T.Trevethan
# 2015
#
# Python script to hydrogenate carbon atomistic structures containing 
# undercoordinated atoms
#
# The structure is supplied as an xyz format file: input.xyz
# 
# The script analyses the input structure and determines the 
# coordination and hybridisation of each carbon atom in the system.
# If a carbon atom has a coordination number of less than 3, a hydrogen
# atom is added to the configuration in a position corresponding to the
# dangling bond (i.e. maximimising the angles to the existing bonds). 
#
# The final hydrogenated configuration is written to the xyz format file:
# output.xyz
#
# The script assumes a C-H bond length of 1.0 A, and a C-C maximum bond
# length of 1.7 A. 
#
# These defaults can be over-ridden with the optional command line arguments:  
#                               
#       -b float : set the C-H bond length
#       -c float : set the C-C bond length
#
# Additionally, optional periodic boundary conditions can be set
#
#       -x float : peridic boundary conditions in x-direction
#       -y float : peridic boundary conditions in y-direction
#       -z float : peridic boundary conditions in z-direction
# 
# Hydrogen atoms placed within 0.7 A of an existing atom are removed
# this tolerance can be over-ridden with the following option:
#
#       -o float : overlap tolerance

import sys
import math 
import argparse

bch = 1.0
bcc = 1.7
overlap = 0.7

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-b', type = float, help = 'C-H bond length')
parser.add_argument('-c', type = float, help = 'C-C bond length')
parser.add_argument('-x', type = float, help = 'x-direction PBC')
parser.add_argument('-y', type = float, help = 'y-direction PBC')
parser.add_argument('-z', type = float, help = 'z-direction PBC')
parser.add_argument('-o', type = float, help = 'overlap tolerance')
variables = parser.parse_args()
if variables.b:
    bch = variables.b
if variables.c:
    bcc = variables.c
if variables.x:
    xper = variables.x
if variables.y:
    yper = variables.y
if variables.z:
    zper = variables.z
if variables.o:
    overlap = variables.o

print "Hydro: program to hydrogenate undercoordinated carbon systems"

# Open input files
input  = open("input.xyz", 'r')

#read input file and store carbon atom coordinates in array
coords = []
label = []

numline = input.readline()
nat = int(numline)
blank = input.readline()

# read coordinates
nc = 0
for i in range(nat):
    line = input.readline()
    fdata = line.split()
    ic = 0
    if fdata[0] == "C":
        ic = 1
        nc += 1
    data = (int(ic),float(fdata[1]),float(fdata[2]),float(fdata[3]))
    coords.append(data)
    label.append(fdata[0])

print "Read in "+str(nat)+" atoms"+" of which "+str(nc)+" are carbon"

#create bond table
bond = []

#added H atoms
hat = []

for i in range(nat):
    ib = 0
    bond.append([])
    for j in range(nat):
        if i != j:
            xdst = coords[i][1] - coords[j][1]
            ydst = coords[i][2] - coords[j][2]
            zdst = coords[i][3] - coords[j][3]
            dst = math.sqrt(xdst*xdst+ydst*ydst+zdst*zdst)
            if dst <= 0.7:
                print "Error: atoms "+str(i)+" and "+str(j)+" closer than 0.7"
                sys.exit(1)
            if dst < bcc:
                if ib > 3:
                    print "Error: atom "+str(i)+" has more than 4 bonds"
                    sys.exit(1)
                bond[i].append(j)
                ib += 1
# if the atom is a carbon
    if coords[i][0] == 1:        
# if the atom is 0 coordinated, leave it alone
# if the atom is 1 coordinated, add an sp H
        if len(bond[i]) == 1:
            xtmp = coords[bond[i][0]][1] - coords[i][1]
            ytmp = coords[bond[i][0]][2] - coords[i][2]
            ztmp = coords[bond[i][0]][3] - coords[i][3]
            bdist = math.sqrt(xtmp*xtmp+ytmp*ytmp+ztmp*ztmp)
            xadd = coords[i][1] - xtmp*bch/bdist
            yadd = coords[i][2] - ytmp*bch/bdist
            zadd = coords[i][3] - ztmp*bch/bdist
            hdata = (float(xadd),float(yadd),float(zadd))
            hat.append(hdata)
# if the atom is 2 coordinated, add an sp2 H
        if len(bond[i]) == 2:
            xt1 = coords[bond[i][0]][1] - coords[i][1]
            yt1 = coords[bond[i][0]][2] - coords[i][2]
            zt1 = coords[bond[i][0]][3] - coords[i][3]
            xt2 = coords[bond[i][1]][1] - coords[i][1]
            yt2 = coords[bond[i][1]][2] - coords[i][2]
            zt2 = coords[bond[i][1]][3] - coords[i][3]
            xadd = coords[i][1] - (xt1+xt2)*bch/1.42
            yadd = coords[i][2] - (yt1+yt2)*bch/1.42
            zadd = coords[i][3] - (zt1+zt2)*bch/1.42
            hdata = (float(xadd),float(yadd),float(zadd))
            hat.append(hdata)

print "Created "+str(len(hat))+" H atoms"
print "Checking for overlap ..."

hatr = 0

# check for overlapping H atoms and remove
for i in range(len(hat)):
    for j in range(len(hat)):
        if j > i:
            xt = hat[i][0] - hat[j][0]
            yt = hat[i][1] - hat[j][1]
            zt = hat[i][2] - hat[j][2]
            dist = math.sqrt(xt*xt+yt*yt+zt*zt)
            if dist < overlap:
                hat[j].append(1)
                hatr += 1

print "Removed "+str(hatr)+" overlapping H atoms"

# open output file and write atoms in xyz format
output = open("output.xyz",'w')

totat = len(coords)+len(hat)-hatr
outline = str(totat)+"\n"
output.write(outline)
output.write(" \n")
for i in range(len(coords)):
    outline = label[i]+" "+str(coords[i][1])+" "+str(coords[i][2])+" "+str(coords[i][3])+" \n"
    output.write(outline)
for i in range(len(hat)):
    if len(hat[i]) == 3:
        outline = "H "+str(hat[i][0])+" "+str(hat[i][1])+" "+str(hat[i][2])+" \n"
        output.write(outline)
