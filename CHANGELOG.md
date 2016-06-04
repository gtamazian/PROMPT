Changes
=======

* 1.2
-----

- Fortran implementation of the transformation cost objective function
added;
- `trmobjfunc` supports planar angles as variables for optimization;
- Functions `trmplotminatomdist` and `pdbplotminatomdist` to plot minimal
interatomic distances; 
- Function `trmchangeangles` for convenient angle change in a model;
- Function `trminteroptim` implementing a multistage optimization
procedure;
- Function `trminitialpoint` to get an initial point for the specified
transformation model and planar and torsion angle indices;
- Functions `pdbstrip` and `pdbenrich` to process PDB files;
- Function `trmanimcost` to produce contour plots of the transformation
  cost function;
- Function `pdbanimrmch` to produce a series of Ramachandran plots for a
  multimodel PDB file;
- Function `pdbcartinterp` to perform linear interpolation of Cartesian
  atom coordinates between two PDB structures;
- Option to leave only alpha carbon atoms added to `pdbbackbone`;
- Option to consider planar angles added to `trmdistantangleindices`;
- `trmcreate` creates an alpha carbon-based transformation model if the
specified PDB structures contains only alpha carbon atoms;
- `trmcreate` checks if the specified PDB structures contain multiple
  chains;
- `pdbcost` supports PDB files containing alpha carbon-based structures;
- `trmplottranglediff` replaced with the `trmplotanglediff` function that
supports both planar and torsion angles;
- Cell arrays of configuration coordinates were replaced with 3D-arrays
for better performance.

1.1
---

- Function `trmcostangles` removed; 
- Plot functions use different markers for lines;
- Y label fixed in plots: 'RMSD in Angstroms' instead of
'RMSD in AA'.

