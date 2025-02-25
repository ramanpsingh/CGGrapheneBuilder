#!/usr/bin/python

import sys
import math
from math import *
import argparse

# Argument Parser
parser = argparse.ArgumentParser(description='Create a custom coordinate and topology file for coarse-grained graphene')
parser.add_argument( "-nt", "--numtriangles", type=int, default=5, help='Number of triangles per row (Corresponds to width of graphene sheet)')
parser.add_argument( "-nr", "--numrows", type=int, default=5, help='Number of rows of triangles (Corresponds to length of graphene sheet)')
parser.add_argument( "-nl", "--numlayers", type=int, default=5, help='Number of layers of graphene sheets (Corresponds to thickness of graphene sheet)')
parser.add_argument( "-bl", "--bondlength", type=float, default=0.47, help='Bond length')
parser.add_argument( "-ang", "--bondangle", type=int, default=60, help='Bond angle')
parser.add_argument( "-bead", "--beadtype", type=str, default='CNP', help='Martini bead type')
parser.add_argument( "-o", "--output",   type=str,   default=None,   help='Name of output file' )
args = parser.parse_args()

lx = args.numtriangles
ly = args.numrows
nl = args.numlayers
l = args.bondlength
ang = args.bondangle
bd = args.beadtype


if args.output == None:
	output = "Graphene-w"+str(lx)+"-l"+str(ly)
else:
        output = args.output

# Number of atoms/beads per row of triangles
natoms_lx = 2*lx+1

# Total number of atoms in ly rows  
natoms_perlayer = natoms_lx*ly

# Total number of atoms
numatoms = natoms_perlayer*nl

# CG graphene is composed of equilateral triangles hence, angle between sides 60 degress
ang=60

# Distance between graphene layers
ldist=0.33

print( "Generated Martini model for Graphene with %s rows of triangles with %s triangles per row " % (ly, lx))

# CREATE GRO FILE

# Opens the gro file for writing

structure_file = open(output+".gro", 'w')


# Writes header for gro file

structure_file.write( "Graphene-w%s-l%s\n" % (lx, ly) )
structure_file.write( "  %3d\n" % (numatoms) )


# Writes coordinates for beads

y_offset = sqrt(3)*l
for k in range(0,nl):
	z = ldist*k
	for j in range(0,ly):
		for i in range(1,natoms_lx+1):
			n= i + natoms_lx*j + natoms_perlayer*k
			x=(i-1)*(l/2)
			if i%2==1:
				y=0 + j*y_offset
			else:
				y=(y_offset/2) + j*y_offset
			
			structure_file.write(  "%5d%-5s G%03d%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "GRA", n, n, x, y, z, 0, 0, 0) )	

# Box dimensions

xdim = lx*l + 1
ydim = y_offset*ly
zdim = ldist*nl + 1

print( "Generated Martini model for Graphene with size %s x %s nm " % (lx*l, y_offset*ly))

structure_file.write(  "   %u  %u  %u\n" % (xdim, ydim, zdim) )

structure_file.close()


# CREATE TOPOLOGY FILE

# Opens the itp file for writing

topology_file = open(output+".itp", 'w')

# Writes header for the itp file

topology_file.write( "; \n; Martini topology for Graphene with with %s rows of triangles with %s triangles per row\n" % (ly, lx))
topology_file.write( "; \n; Topology generated using GrapheneCGBuilder.py \n" )
topology_file.write( "; Written by Raman Preet Singh \n\n")

topology_file.write( "[ moleculetype ]\n" )
topology_file.write( "; Name	 nrexcl\n" )
topology_file.write( "GRA  1\n")

# Atoms

topology_file.write( "\n[ atoms ]\n" )
topology_file.write( "; id	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n" )

for k in range(0,nl):
	for j in range(0,ly):
		for i in range(1,natoms_lx+1):
			n= i + natoms_lx*j + natoms_perlayer*k
			topology_file.write( "%3d    %4s   1   GRA    G%03d     %3d     0       48\n" % (n, bd, n, n) )

# Bonds

topology_file.write( "\n[ bonds ]\n" )
topology_file.write( "; i	 j	  funct	 length	 force\n" )

## Bonds within a layer

for k in range(0,nl):
	for j in range(0,ly):
		for i in range(1,natoms_lx-1):
			m = i + natoms_lx*j + natoms_perlayer*k
			n = m + 1
			o = m+2
			topology_file.write( "     %3d     %3d       1   %4.3f     10000\n" % (m, n, l) )
			topology_file.write( "     %3d     %3d       1   %4.3f     10000\n" % (m, o, l) )
		else:	
			m = i + natoms_lx*j + natoms_perlayer*k + 1
			n = m + 1
			topology_file.write( "     %3d     %3d       1   %4.3f     10000\n" % (m, n, l) )

for k in range(0,nl):
	for j in range(0,ly-1):
		for i in range(1,natoms_lx):
			if i%2==0:
				m = i + natoms_lx*j + natoms_perlayer*k
				n = i + natoms_lx*(j+1) - 1 + natoms_perlayer*k
				o = n + 2
				topology_file.write( "     %3d     %3d       1   %4.3f     25000\n" % (m, n, l) )
				topology_file.write( "     %3d     %3d       1   %4.3f     25000\n" % (m, o, l) )

## Bonds between layers
for k in range(0,nl-1):
	for j in range(0,ly):
		for i in range(1,natoms_lx+1):
			m= i + natoms_lx*j + natoms_perlayer*k
			n = i + natoms_lx*j + natoms_perlayer*(k+1)
			topology_file.write( "     %3d     %3d       1   %4.3f     25000\n" % (m, n, ldist) )


# Angles

topology_file.write( "\n[ angles ]\n" )
topology_file.write( "; i	 j	 k	 funct	 angle	 force\n" )


## Angles within a layer
for k in range(0,nl):
	for j in range(0,ly):
		for i in range(1,natoms_lx-1):
			m = i + natoms_lx*j + natoms_perlayer*k
			n = m + 1
			o = m+2
			p = m+4
			q = m + natoms_lx
			r = q+2
			s=q+1
			topology_file.write( "     %3d     %3d     %3d       2     %3.1f     500.0\n" % (m, n, o, ang) )
			topology_file.write( "     %3d     %3d     %3d       2     %3.1f     500.0\n" % (n, m, o, ang) )
			topology_file.write( "     %3d     %3d     %3d       2     %3.1f     500.0\n" % (n, o, m, ang) )
			if p<=natoms_lx*(j+1):
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (m, o, p, 180) )
			if i%2!=0 and r<=natoms_perlayer*(k+1):
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (m, n, q, 2*ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (o, n, r, 2*ang) )
			if i%2==0 and r<=natoms_perlayer*(k+1):
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (m, s, q, 2*ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (o, s, r, 2*ang) )

for k in range(0,nl):
	for j in range(0,ly-1):
		for i in range(1,natoms_lx):
			if i%2==0:
				m = i + natoms_lx*j + natoms_perlayer*k
				n = i + natoms_lx*(j+1) - 1
				o = n + 2
				p = m + 2
				q = m + 1
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (m, n, o, ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (n, m, o, ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (n, o, m, ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (o, m, p, ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (o, p, m, ang) )
				topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (m, o, p, ang) )

## Angles between layers
for k in range(0,nl-1):
	for j in range(0,ly):
		for i in range(1,natoms_lx-1):
			m= i + natoms_lx*j + natoms_perlayer*k
			n = m + 1
			o = i + natoms_lx*j + natoms_perlayer*(k+1) + 1
			p = i + natoms_lx*j + natoms_perlayer*(k+1)
			topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (m, n, o, 90) )
			topology_file.write( "     %3d     %3d     %3d       2     %3.1f    500.0\n" % (n, m, p, 90) )


# Dihedrals

topology_file.write( "\n[ dihedrals ]\n" )
topology_file.write( "; i	 j	 k	 l     func	 q0     cq\n" )

for k in range(0,nl):
	for j in range(0,ly-1):
		for i in range(2,natoms_lx-1):
			m = i + natoms_lx*j + natoms_perlayer*k
			n = m + 1
			o = m + 3
			p = m + 4
			q = m + natoms_lx + 1
			r = m + natoms_lx + 3
			s = m + natoms_lx
			if i%2==0:
				topology_file.write( "     %3d     %3d     %3d     %3d       9     %3.1f     500.0\n" % (n, m, q, r, 180) )
				if p <= (natoms_lx*(j+1) + natoms_perlayer*k):
					topology_file.write( "     %3d     %3d     %3d     %3d       9     %3.1f     500.0\n" % (m, n, o, p, 180) )
			else:
				
				if p <= (natoms_lx*(j+1) + natoms_perlayer*k):
					topology_file.write( "     %3d     %3d     %3d     %3d       9     %3.1f     500.0\n" % (n, s, q, r, 180) )



# Improper dihedrals
topology_file.write( "\n[ dihedrals ]\n" )
topology_file.write( "; i	 j	 k	 l     func	 q0     cq\n" )

for k in range(0,nl):
	for j in range(0,ly-1):
		for i in range(2,natoms_lx-1):
			if i%2==0:
				m = i + natoms_lx*j + natoms_perlayer*k
				n = m + 1
				o = m+2
				p = m +1 + natoms_lx
				topology_file.write( "     %3d     %3d     %3d     %3d       2     %3.1f     500.0\n" % (n, m, o, p, 180) )


for k in range(0,nl):
	for j in range(0,ly-1):
		for i in range(2,natoms_lx):
			if i%2==0:
				m = i + natoms_lx*j + natoms_perlayer*k
				n = m + natoms_lx - 1
				o = m + natoms_lx + 1
				p = m + natoms_lx
				topology_file.write( "     %3d     %3d     %3d     %3d       2     %3.1f     500.0\n" % (m, n, o, p, 180) )

# Restraints

topology_file.write( "\n; Include Position restraint file\n" )
topology_file.write( "#ifdef POSRES\n" )
topology_file.write( "#include \"%s-posres.itp\"\n" % (output) )
topology_file.write( "#endif\n" )

topology_file.close()



# CREATE POSITION --------------------------#
# Position Restraints File #
#--------------------------#


# Open the file for writing

posres_file = open(output+"-posres.itp", 'w')


# Header

posres_file.write( "[ position_restraints ]\n" )
posres_file.write( "; ai  funct  fcx    fcy    fcz\n" )

for i in range(1, numatoms+1):
	posres_file.write( " %3d    1    1000   1000   1000\n" % (i) )
	
posres_file.close()

