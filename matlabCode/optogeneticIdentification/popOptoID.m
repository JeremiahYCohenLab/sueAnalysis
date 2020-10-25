function [spikeLat, spikeProb] = popOptoID(xlFile, sheet)

[~, sessionCellList, ~] = xlsread(xlFile, sheet);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);

for i = 1:length(cellList)
    sessionList(i)
    cellList(i)
    spikeLatTmp = optoIDnlx(sessionList{i}, 'cellName', cellList{i}, 'plotFlag', 0);
    spikeLatTmp = spikeLatTmp(:);
    if spikeLatTmp == 0
        spikeLat(i) = NaN;
        spikeProb(i) = NaN;     
    elseif any(~isnan(spikeLatTmp))
        spikeLat(i) = nanmean(spikeLatTmp);
        spikeProb(i) = sum(sum(~isnan(spikeLatTmp)))/100;
    else
        spikeLat(i) = 0;
        spikeProb(i) = 0;
    end
end
spikeLat(isnan(spikeLat)) = [];
spikeProb(isnan(spikeProb)) = [];
figure; hold on;
scatter(spikeLat, spikeProb);

end