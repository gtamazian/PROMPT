% PROMPT Toolbox
% Version 1.0 15-Dec-2014
%
% PROMPT (PRotein cOnformation Motion PredicTion) toolbox contains a set 
% of functions to create and process coarse-grained models of protein 
% conformational motion. The method implemented in the toolbox is described
% in the following paper.
%
% Modeling of Conformation Change by Redox-Switch Modulation of Human
% Succinic Semialdehyde Dehydrogenase by Gaik Tamazian, Jeong Ho Chang, 
% Sergey Knyazev, Eugene Stepanov, Kyung-Jin Kim and Yuri Porozov.
% Submitted to PLoS Computational Biology journal.
%
% Demos
%   calmodulin_demo - Calmodulin conformational motion simulation demo
%
% Functions
%   atomicmass               - Get atomic masses
%   circdist                 - Calculate circular distance between angles
%   circinterp               - Circular interpolation of angles by the specified arcs
%   createmodel              - Get transformation bond lengths, planar and torsion angles
%   cross3d                  - Simplified evalulation of a vector cross product for 3D vectors
%   getbinarycode            - Generate binary numbers from 0 to 2^N
%   optimquat                - Get quaternion representing the optimal rotation to fit Y to X
%   pdb2aa                   - Return a cell array of PDB structure amino acid sequences
%   pdb2trm                  - Create a transformation model from a PDB structure
%   pdbbackbone              - Return a PDB structure only with backbone atoms
%   pdbextractcoords         - Extract Cartesian atom coordinates from a PDB structure
%   pdbmininteratomicdist    - Returns minimal interatomic distances
%   pdbplotadjrmsd           - Plot RMSDs between adjacent models
%   pdbplotfixedrmsd         - Plot RMSDs between models and the specified one
%   pdbtrfcmp                - Compare transformations from PDB structures
%   pdbtrfcost               - Calculate transformation cost from a PDB structure
%   quat2rotmat              - Convert a quaternion to a rotation matrix
%   reduceangles             - Reduce angular values to the interval [-pi, pi]
%   restorecoords            - Restore Cartesian coordinates of atoms
%   sidechainmass            - Get atomic masses of protein side chains
%   sievefit                 - Apply the sieve-fit procedure to fit Y to X
%   trm2pdb                  - Convert a transformation model to a PDB structure object
%   trmcost                  - Calculate transformation cost
%   trmcostangles            - Transformation cost as a function of specified angles
%   trmcreate                - Create a transformation model
%   trmdistantangleindices   - Get indices of the most distant torsion angles
%   trmmininteratomicdist    - Get minimal interatomic distances
%   trmobjfunc               - Calculate the cost function and its gradient for the model
%   trmplotadjrmsd           - Plot RMSDs between adjacent configurations
%   trmplotfixedrmsd         - Plot RMSDs between configurations and the specified one
%   trmplottranglediff       - Plot differences between transformation torsion angles
%   trmrestorecoords         - Restore Cartesian coordinates of transformation atoms
%   trmupdaterotations       - Update rotation matrices by the optimal superposition
%   trmvariateinterpolations - Vary interpolation modes for torsion angles
