function n = initModelMaxInd(firstConfPdbPath, lastConfPdbPath, ...
    firstConfInd, lastConfInd, bidirAngThrhld)

bidirAngThrhld = double(bidirAngThrhld);

firstConf = pdbread(firstConfPdbPath);
firstConf.Model = firstConf.Model(firstConfInd);

lastConf = pdbread(lastConfPdbPath);
lastConf.Model = lastConf.Model(lastConfInd);

trm = trmcreate(firstConf, lastConf, 0);
n = 2^sum(abs(circdist(trm.psi(:,1), trm.psi(:,end))) >= bidirAngThrhld);
end
