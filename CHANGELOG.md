Changes
=======

* 1.2
-----

- Fortran implementation of the transformation cost objective function
added;
- Option to leave only alpha carbon atoms added to `pdbbackbone`;
- `trmcreate` creates a alpha carbon-based transformation model if the
specified PDB structures contains only alpha carbon atoms;
- `trmplottranglediff` replaced with the `trmplotanglediff` function that
supports both planar and torsion angles.

1.1
---

- Function `trmcostangles` removed; 
- Plot functions use different markers for lines;
- Y label fixed in plots: 'RMSD in Angstroms' instead of
'RMSD in AA'.

