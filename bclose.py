#!/usr/bin/python

# bclose.py
# T.Trevethan
# 2016
#
# Python script perform a bond distance analysis on a supllied atomistic
# configuration. The code produces a radial distribution function and 
# can additionally remove atoms that are closer than a supplied cutoff. 
# 
# The link cell algorithm is employed to determine neighbour lists, which
# means that even very large systems (millions of atoms) can be analysed
# quickly and efficiently
#
# The structure is supplied as an xyz format file: input.xyz
# 
# If the remove close option is selected, the final corrected 
# configuration is written to the xyz format file:output.xyz. This is set
# (along with the cutoff) using:
#
#       -c float: remove atoms within the supplied cutoff
#
# The RDF is writen to the file rdf.dat. The default precision is 0.1 A 
# and the range is 10 A. These can be over-ridden with the optional
# command line arguments:
#
#       -r float: set the RDF resolution
#       -m float: set the maxium of the RDF
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
# The link-cell devisor is set by default to 16. To over-ride this set:
#
#       -i int: link cell devisor

import sys
import math 
import argparse

# function determine distance between two coordinates
def dist(coords1,coords2):
    x1 = coords1[0] - coords2[0]
    y1 = coords1[1] - coords2[1]
    z1 = coords1[2] - coords2[2]
    bdist = math.sqrt(x1*x1+y1*y1+z1*z1)
    return bdist

res = 0.1
rmax = 10.0
icdim = 16

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-r', type = float, help = 'RDF resoltuion')
parser.add_argument('-m', type = float, help = 'RDF maximum')
parser.add_argument('-x', type = float, help = 'x-direction PBC')
parser.add_argument('-y', type = float, help = 'y-direction PBC')
parser.add_argument('-z', type = float, help = 'z-direction PBC')
parser.add_argument('-c', type = float, help = 'remove close tolerance')
parser.add_argument('-i', type = int, help = 'link cell divisor')
variables = parser.parse_args()
if variables.r:
    res = variables.r
if variables.m:
    rmax = variables.m
if variables.x:
    xper = variables.x
if variables.y:
    yper = variables.y
if variables.z:
    zper = variables.z
if variables.c:
    rclose = variables.c
if variables.i:
    icdim = variables.i

# initialise histogram bins
rdf = []
nbins = int(rmax/res)
for i in range(nbins):
    rdf.append(0)

print "Bclose: program to perform a bond neighbour analysis and remove close atoms"

# Open input file
input  = open("input.xyz", 'r')

#read input file and store atom coordinates in array
coords = []
label = []

numline = input.readline()
nat = int(numline)
blank = input.readline()

# read coordinates
# and determine the coordinate ranges
xmin = ymin = zmin = 0.0
xmax = ymax = zmax = 0.0
for i in range(nat):
    line = input.readline()
    fdata = line.split()
    if float(fdata[1]) < xmin:
        xmin = float(fdata[1])
    elif float(fdata[1]) > xmax:
        xmax = float(fdata[1])
    if float(fdata[2]) < ymin:
        ymin = float(fdata[2])
    elif float(fdata[2]) > ymax:
        ymax = float(fdata[2])
    if float(fdata[3]) < zmin:
        zmin = float(fdata[3])
    elif float(fdata[3]) > zmax:
        zmax = float(fdata[3])
    data = (float(fdata[1]),float(fdata[2]),float(fdata[3]))
    coords.append(data)
    label.append(fdata[0])

print "Read in "+str(nat)+" atoms"

#get system dimensions
if xper:
    print "PBC in x direction"
    xln = xper
else:
    print "Min x: "+str(xmin)+" Max x: "+str(xmax)
    xln = xmax - xmin

if yper:
    print "PBC in y direction"
    yln = yper
else:
    print "Min y: "+str(ymin)+" Max y: "+str(ymax)
    yln = ymax - ymin

if zper:
    print "PBC in z direction"
    zln = zper
else:
    print "Min z: "+str(zmin)+" Max z: "+str(zmax)
    zln = zmax - zmin

icell = []
jcell = []
kcell = []
# create link cell lists
for k in range(icdim):
    kcell.append(0)

for j in range(icdim):
    jcell.append(kcell)

for i in range(icdim):
    icell.append(jcell)

xdiv = xln/(icdim*1.0)
ydiv = xln/(icdim*1.0)
zdiv = xln/(icdim*1.0)

if xdiv < rmax or ydiv < rmax or zdiv < rmax:
    print "Error: RDF max greater than link cell size"
    print "Either decrease RDF max or cell divisor"
    sys.exit(1)

for i in xrange(len(coords)):
    icx = int(coords[i][0]/xdiv)
    icy = int(coords[i][1]/ydiv)
    icz = int(coords[i][2]/zdiv)
    if icx > icdim:
        icx = icdim
    if icy > icdim:
        icy = icdim
    if icz > icdim:
        icz = icdim
    if icx < 1:
        icx = 1
    if icy < 1:
        icy = 1
    if icz < 1:
        icz = 1

    icell[icx-1][icy-1][icz-1].append(i)

# removal list
rlist = []

tcrd = []
for i in range(3):
    tcrd.append(0.0)

# loop over all cells to evaluate neighbour distances
for i in range(icdim):
    for j in range(icdim):
        for k in range(icdim):
#           double loop over all atoms in the cell
            for n in range(len(icell[i][j][k])):
                iatm1 = icell[i][j][k][n]
                for m in range(len(icell[i][j][k])):
                    if m != n:
                        iatm2 = icell[i][j][k][m]
                        bdst = dist(coords[iatm1],coords[iatm2])
                        # add to remove list if too close to neighbour
                        if rclose:
                            if bdst < rclose:
                                rlist.append(iatm2)
                        # add neoighbour distance to RDF histogram
                        for p in xrange(nbins):
                            if bdst > p*res and bdst < (p+1)*res:
                                rdf[p] += 1
#    loop over all neighbouring cells
                for icp in range(-1,1):
                    for jcp in range(-1,1):
                        for kcp in range(-1,1):
                            if icp==0 and jcp==0 and kcp==0:
                                break
#               do periodic boundary conditions
                            xaddp = 0.0
                            yaddp = 0.0
                            zaddp = 0.0
                            iper = i + icp
                            jper = j + jcp
                            kper = k + kcp
                            if iper == -1:
                                if xper:
                                    iper = icdim
                                    xaddp = -xln
                                else:
                                    break
                            if jper == -1:
                                if yper:
                                    jper = icdim
                                    yaddp = -yln
                                else:
                                    break
                            if kper == -1:
                                if zper:
                                    kper = icdim
                                    zaddp = -zln
                                else:
                                    break
                            if iper = icdim:
                                if xper:
                                    iper = 0
                                    xaddp = xln
                                else:
                                    break
                            if jper = icdim:
                                if yper:
                                    jper = 0
                                    yaddp = yln
                                else:
                                    break
                            if kper = icdim:
                                if zper:
                                    kper = 0
                                    zaddp = zln
                                else:
                                    break
#         loop over atoms in cell
                            for m in range(icell[iper][jper][kper]):
                                iatm2 = icell[iper][jper][kper][m]
                                tcrd[0] = coords[iatm2][0] + xaddp
                                tcrd[1] = coords[iatm2][1] + yaddp
                                tcrd[2] = coords[iatm2][2] + zaddp
                                bdst = dist(coords[iatm1],tcrd)
                                if rclose:
                                    if bdst < rclose:
                                        rlist.append(iatm2)
# add neoighbour distance to RDF histogram
                                for p in xrange(nbins):
                                    if bdst > p*res and bdst < (p+1)*res:
                                        rdf[p] += 1

if rclose:
    print "Removed "+str(len(rlist))+" close atoms"

#  write coordinates to output minus the removed atoms
    output = open("output.xyz",'w')

    totat = len(coords)-len(rlist)
    outline = str(totat)+"\n"
    output.write(outline)
    output.write(" \n")
    for i in range(len(coords)):
        skip = 0
        for j in range(len(rlist)):
            if rlist[j] == i:
                skip = 1
        if skip == 0:
            outline = label[i]+" "+str(coords[i][0])+" "+str(coords[i][1])+" "+str(coords[i][2])+" \n"
            output.write(outline)

# write RDF histogram
outrdf = open("rdf.dat",'w')

for i in range(nbins):
    rdfline = str(res*i)+" "+str(rdf[i])
    outrdf.write(rdfline)
