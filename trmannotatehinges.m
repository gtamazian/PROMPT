function trm = trmannotatehinges(trm, minPhi, maxGap, maxGapToTail)
%TRMANNOTATEHINGES Finds and adds an annotation about hinges to trm

% By Sergey Knyazev, 2016.
% sergey (dot) n (dot) knyazev (at) gmail (dot) com

MIN_PHI = 0.5;
MAX_GAP = 30;
MAX_GAP_TO_TAIL = 30;

if ~exist('minPhi', 'var') || ~minPhi
    minPhi = MIN_PHI;
end

if ~exist('maxGap', 'var') || ~maxGap
    maxGap = MAX_GAP;
end

if ~exist('maxGapToTail', 'var') || ~maxGapToTail
    maxGapToTail = MAX_GAP_TO_TAIL;
end

angdiff = abs(circdist(trm.psi(:,1), trm.psi(:,end)));
ind = 1:length(angdiff);
hingeAngles = ind(angdiff > minPhi);

if ~isempty(hingeAngles)
    gaps = hingeAngles(2:end) - hingeAngles(1:end-1);
    ind = 1:length(gaps);
    gapsIndices = ind(gaps > maxGap);
    trm.hinges = [hingeAngles([1 (gapsIndices + 1)]).' ...
        hingeAngles([gapsIndices length(hingeAngles)]).'];
    trm.hinges = trm.hinges(trm.hinges(:,1) ~= trm.hinges(:,2),:);

    trm.hingeTypes = repmat({'hinge'}, size(trm.hinges,1), 1);
    if trm.hinges(1,1) <= maxGapToTail
        trm.hingeTypes{1} = 'n-tail';
        trm.hinges(1,1) = 1;
    end
    if trm.hinges(end,end) >= size(trm.psi,1) - maxGapToTail
        trm.hingeTypes{end} = 'c-tail';
        trm.hinges(end,end) = size(trm.psi,1);
    end
end

end

