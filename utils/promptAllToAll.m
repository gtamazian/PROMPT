function promptAllToAll(inputPdbPath, outPdbPath, ...
    nConf, nAngles, bidirAngThrhld, iterNum, superpose)

pdbStructure = pdbread(inputPdbPath);
firstConf = pdbStructure;
lastConf = pdbStructure;
n = length(pdbStructure.Model);
for i=1:n
    firstConf.Model = pdbStructure.Model(i);
    for j = i+1:n
        lastConf.Model = pdbStructure.Model(j);
        m = initModelMaxInd(firstConf, lastConf, bidirAngThrhld);
        for k = 1:m
        model = prompt(firstConf, lastConf, nConf, nAngles, ...
            k, iterNum, superpose);
        pdbOptimized = trm2pdb(model, pdbStructure);
        pdbwrite([outPdbPath int2str(i) '_' int2str(j) '_' ...
            int2str(k) '.pdb'], ...
            pdbOptimized);
        end
    end
end
end

function n = initModelMaxInd(firstConf, lastConf, bidirAngThrhld)
trm = trmcreate(firstConf, lastConf, 0);
n = 2^sum(abs(circdist(trm.psi(:,1), trm.psi(:,end))) >= bidirAngThrhld);
end
