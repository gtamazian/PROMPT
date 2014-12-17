%% Calmodulin conformational motion simulation demo
% We demonstrate routines of the PROMPT package by modeling
% conformational motion between two conformations of calmodulin. To launch
% the example, use the following command: *calmodulin_demo*.

%% Preparing data
% Prior to creating a transformation model, we should process data from PDB
% files. First, we load data from the original calmodulin PDB file using
% the *pdbread* function from _MATLAB Bioinformatics Toolbox_
% (<http://www.mathworks.com/products/bioinfo/>). The function returns a
% *PDBStruct* object.

toolboxPath = fileparts(which('calmodulin_demo'));
pdbStructure = pdbread(fullfile(toolboxPath, 'samples/calmodulin.pdb'));

%%
% We need a pair of protein states to build a transformation model. Further
% we call any protein state a _configuration_. The configurations are 
% exracted from the obtained PDBStruct object in the following way.

firstConf = pdbStructure; firstConf.Model = firstConf.Model(1);
lastConf = pdbStructure; lastConf.Model = lastConf.Model(end);

%%
% The obtained PDBStruct objects *firstConf* and *lastConf* contain a 
% single configuration and will be used to create a transformation model.

%% Creating a transformation model
% We create a transformation model from a pair of *PDBStruct* objects using
% the *trmcreate* function from the _PROMPT_ toolbox. The number of
% intermediate configurations is also specified as an argument of the
% function.

nConf = 8;
model = trmcreate(firstConf, lastConf, nConf);

%%
% The next step is to optimize the created model. But prior to the
% optimization, we study the model in order to determine appropriate
% optimization parameters. First, we visualise the difference between the 
% first and last configurations using the *trmplottranglediff* function.

trmplottranglediff(model, true);
grid;
title(['Difference Between Torsion Angles of the First and ', ...
    'Last Configurations']);

%%
% Less than 50 angles differ in more than 20 degress. Let us plot the
% differences for 50 angles which differ the most.

nAngles = 50;
trmplottranglediff(model, true, nAngles, '-o');
grid;
title(['50 Torsion Angles Differing the Most between the First and ', ...
    'Last Configurations']);

%%
% Let us choose 20 torsion angles for the optimization procedure.

nAngles = 20;

%% Optimization of the transformation model
% We optimize the transformation model by the interior-point algorithm
% implemented in the *fmincon* function from _MATLAB Optimization Toolbox_.
% First, we create the objective function *f* to be passed to *fmincon*
% using the *trmobjfunc* function from the _PROMPT_ package.

angleIndices = trmdistantangleindices(model, nAngles);
f = @(x) trmobjfunc(model, angleIndices, x);

%%
% Next, we form the initial point to start optimization from.

initial_point = model.psi(angleIndices,2:end-1);
initial_point = initial_point(:);
initial_point = reduceangles(initial_point);

%%
% Also we specify the optimization parameters. The parameter *iterNum*
% defines the maximum number of optimization procedure iterations. In this
% demonstration, we specify its value equal to 100 for fast result 
% computation.

iterNum = 100;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxIter', iterNum);
options = optimoptions(options,'MaxFunEvals', Inf);
options = optimoptions(options,'GradObj','on');

%%
% Finally, we launch the optimization procedure.

n = size(initial_point);
x = fmincon(f,initial_point,[],[],[],[],-pi*ones(n),pi*ones(n),[],...
    options);

%%
% To store the optimized transformation, we create a separate model.
modelOptimized = model;
modelOptimized.psi(angleIndices, 2:end-1) = reshape(x, ...
    length(angleIndices), size(modelOptimized.psi, 2) - 2);

%% Studying the optimization result
% To compare the optimized transformation with the original one, we plot
% RMSDs between its configurations. RMSDs between adjacent configurations
% are plotted using the function *trmplotadjrmsd*.

trmplotadjrmsd({model, modelOptimized});
legend('original', 'optimized');
grid;
title('RMSDs Between Adjacent Configurations');

%%
% RMSDs to the fixed configuration (we choose the first one) are plotted 
% using the function *trmplotfixedrmsd*.

confNo = 1;
trmplotfixedrmsd({model, modelOptimized}, confNo);
legend('original', 'optimized');
grid;
title('RMSDs to the First Configuration');

%%
% Both RMSDs between the adjacent contigurations and to the first
% configuration decreased after optimization. The larger number of
% optimization itetations *iterNum* will lead to less RMSD values and
% smoother protein motion.

