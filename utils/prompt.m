function cost = prompt(firstConfPdbPath, lastConfPdbPath, ...
    firstConfInd, lastConfInd, outPdbPath, ...
    nConf, nAngles, initPointInd, iterNum, superpose)

firstConf = pdbread(firstConfPdbPath);
firstConf.Model = firstConf.Model(firstConfInd);
firstConf = pdbbackbone(firstConf);

lastConf = pdbread(lastConfPdbPath);
lastConf.Model = lastConf.Model(lastConfInd);
lastConf = pdbbackbone(lastConf);

model = trmcreate(firstConf, lastConf, nConf, initPointInd);

angleIndices = trmdistantangleindices(model, nAngles);
f = @(x) trmobjfunc(model, angleIndices, x, superpose);

initial_point = model.psi(angleIndices,2:end-1);
initial_point = initial_point(:);
initial_point = reduceangles(initial_point);

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxIter', iterNum);
options = optimoptions(options,'MaxFunEvals', Inf);
options = optimoptions(options,'GradObj','on');

n = size(initial_point);
x = fmincon(f,initial_point,[],[],[],[],-pi*ones(n),pi*ones(n),[],...
    options);

modelOptimized = model;
modelOptimized.psi(angleIndices, 2:end-1) = reshape(x, ...
    length(angleIndices), size(modelOptimized.psi, 2) - 2);
if superpose
    modelOptimized.U = trmsuperpos(modelOptimized);
end

cost = trmcost(modelOptimized);
pdb = trm2pdb(modelOptimized, firstConf);
pdbwrite(outPdbPath, pdb);

end

