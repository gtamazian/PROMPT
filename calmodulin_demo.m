%% Calmodulin conformational motion simulation demo
% We demonstrate routines of the _PROMPT_ toolbox by modeling
% conformational motion between two calmodulin structures. To launch
% the example, use the following command: *calmodulin_demo*.

%% Preparing data
% Before creating the transformation model, we should process data from PDB
% files. First, we load data from the original calmodulin PDB file using
% the *pdbread* function from _MATLAB Bioinformatics Toolbox_
% (<http://www.mathworks.com/products/bioinfo/>). The function returns a
% *PDBStruct* object.

toolboxPath = fileparts(which('calmodulin_demo'));
pdbStructure = pdbread(fullfile(toolboxPath, 'samples/calmodulin.pdb'));

%%
% We need a pair of protein structures to build the transformation model. Further
% we call any protein structure a _configuration_. The configurations are
% exracted from the obtained PDBStruct object in the following way.

firstConf = pdbStructure; firstConf.Model = firstConf.Model(1);
lastConf = pdbStructure; lastConf.Model = lastConf.Model(end);

%%
% The obtained PDBStruct objects *firstConf* and *lastConf* contain a 
% single configuration and will be used to create the transformation model.

%% Creating the transformation model
% We create the transformation model from a pair of *PDBStruct* objects using
% the *trmcreate* function from the _PROMPT_ toolbox. The number of
% intermediate configurations is also specified as an argument of the
% function.

nConf = 8;
model = trmcreate(firstConf, lastConf, nConf);

%%
% The next step is to optimize the created model. For that purpose,
% we study the model in order to determine appropriate
% optimization parameters. First, we visualise the difference between the 
% first and last configurations using the *trmplottranglediff* function
% from the _PROMPT_ toolbox.

trmplotanglediff(model, 'torsion', true);
grid;
title(['Difference Between Torsion Angles of the First and ', ...
    'Last Configurations']);

%%
% Less than 50 angles differ in more than 20 degress. Let us plot the
% differences for 50 angles which differ the most.

nAngles = 50;
trmplotanglediff(model, 'torsion', true, nAngles, '-o');
grid;
title(['50 Torsion Angles Differing the Most between the First and ', ...
    'Last Configurations']);

%%
% Let us choose 20 torsion angles for the optimization procedure.

nAngles = 20;

%% Optimization of the transformation model
% We optimize the transformation model by the quasi-Newton algorithm
% implemented in the *fminunc* function from _MATLAB Optimization Toolbox_.
% First, we create the objective function *f* to be passed to *fminunc*
% using the *trmobjfunc* function from the _PROMPT_ toolbox.

angleIndices = trmdistantangleindices(model, nAngles);
f = @(x) trmobjfunc(model, [], angleIndices, x);

%%
% Next, we form the initial point to start optimization from.

initial_point = trminitialpoint(model, [], angleIndices);

%%
% Also we specify the optimization parameters. The parameter *iterNum*
% defines the maximum number of optimization procedure iterations. In this
% demonstration, we specify its value equal to 100 for rapid result
% computation.

iterNum = 100;
options = optimoptions('fminunc', ...
    'Algorith', 'quasi-newton', ...
    'Display', 'iter', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxIterations', iterNum, ...
    'MaxFunctionEvaluations', Inf);
%%
% Finally, we launch the optimization procedure.

n = size(initial_point);
x = fminunc(f,initial_point,options);

%%
% To store the optimized transformation, we create a separate model.
modelOptimized = trmchangeangles(model, [], angleIndices, x);

%% Studying the optimization result
% To compare the optimized transformation with the original one, we plot
% root-mean-square deviations (RMSDs) between its configurations. RMSDs
% between adjacent configurations are plotted using the *trmplotadjrmsd*
% function from the _PROMPT_ toolbox.

trmplotadjrmsd({model, modelOptimized});
legend('original', 'optimized', 'Location', 'best');
grid;
title('RMSDs Between Adjacent Configurations');

%%
% RMSDs to the fixed configuration (we choose the first one) are plotted 
% using the *trmplotfixedrmsd* function from the _PROMPT_ toolbox.

confNo = 1;
trmplotfixedrmsd({model, modelOptimized}, confNo);
legend('original', 'optimized', 'Location', 'best');
grid;
title('RMSDs to the First Configuration');

%%
% Both RMSDs between the adjacent contigurations and to the first
% configuration decreased after optimization. The larger number *iterNum*
% of optimization itetations  will lead to lesser RMSD values and
% smoother simulated motion of the protein.

