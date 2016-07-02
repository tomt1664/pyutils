# pyutils
Collection of python scripts for data processing and structure analysis

> hydro.py

Python script to hydrogenate carbon atomistic structures containing undercoordinated atoms

The structure is supplied as an xyz format file: input.xyz

The script analyses the input structure and determines the coordination and hybridisation of each carbon atom in the system.
If a carbon atom has a coordination number of less than 3, a hydrogenatom is added to the configuration in a position corresponding to the dangling bond (i.e. maximimising the angles to the existing bonds). 

The final hydrogenated configuration is written to the xyz format file: output.xyz

> ringo.py

Script to analyse a covalently bonded structure and output all n-fold bonded rings

The structure is supplied as an xyz format file: input.xyz

The script analyses the input structure and prints the number of 3-fold (triangles) 4-fold (squares) 5-fold (pentagons) 6-fold (hexagons) and 7-fold (heptagons)

The coordinates of the polygons are output to xyz format files if specified at the command line. 

Note: requires numpy

> bclose.py

Python script perform a bond distance analysis on a supllied atomistic configuration. The code produces a radial distribution function and can additionally remove atoms that are closer than a supplied cutoff. 

The link cell algorithm is employed to determine neighbour lists, which means that even very large systems (100s of millions of atoms) can be analysed quickly and efficiently

The structure is supplied as an xyz format file: input.xyz

If the remove close option is selected, the final corrected configuration is written to the xyz format file:output.xyz. 
