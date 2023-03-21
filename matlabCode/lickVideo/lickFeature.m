function [lickS, decisionLicks] = lickFeature(session)
% load data
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
savepath = [pd.sortedFolder session '_tongue.mat'];
load(savepath);
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);

% lick peak detection 
thresh = 0; % threshold for crossing lick port
preBins = 4; % num of bins before the peak for a peak speed to be calculated


ind = 1:length(s.allChoices);
lickS.meanX = NaN(1, length(ind));
lickS.meanY = NaN(1, length(ind));
lickS.meanXPre = NaN(1, length(ind));
lickS.XMax = NaN(1, length(ind));
lickS.YMax = NaN(1, length(ind));
lickS.speedAtPort = NaN(1, length(ind));
lickS.vXAtPort = NaN(1, length(ind));
lickS.vYAtPort = NaN(1, length(ind));
lickS.vXPrePeak = NaN(1, length(ind));
lickS.vYPrePeak = NaN(1, length(ind));
lickS.vXFixInter = NaN(1, length(ind));
lickS.vYFixInter = NaN(1, length(ind));
lickS.speedPrePeak = NaN(1, length(ind));
lickS.crossInd = NaN(1, length(ind));
lickS.maxIndAll = NaN(1, length(ind));
lickS.crossX = NaN(1, length(ind));
lickS.crossY = NaN(1, length(ind));
lickS.len = NaN(1, length(ind));
lickS.conf = NaN(1, length(ind));
decisionLicks = cell(length(ind), 1);
for j = 1:length(ind)
    if allLicks(ind(j)).decisionID
        decisionID = allLicks(ind(j)).decisionID;
        interval = allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1);
        disRL = norm(allLicks(ind(j)).portL - allLicks(ind(j)).portR);
        temp = allLicks(ind(j)).allLicksFilt{decisionID};
        tempX = (temp(:,1) - allLicks(ind(j)).midPoint(1)) * 4/disRL;
        tempY = -(temp(:,2) - allLicks(ind(j)).midPoint(2)) * 4/disRL;
        [meanYMax(ind(j)), maxInd] = max(tempY);
        lickS.maxIndAll(ind(j)) = maxInd;
        lickS.XMax(ind(j)) = tempX(maxInd);
        lickS.YMax(ind(j)) = tempY(maxInd);
        lickS.XPre(ind(j)) = mean(tempX(1:maxInd));
        lickS.meanX(ind(j)) = mean(tempX);
        lickS.meanY(ind(j)) = mean(tempY);
        time = allLicks(ind(j)).time(allLicks(ind(j)).windows(decisionID,1): allLicks(ind(j)).windows(decisionID,2));
        decisionLicks{ind(j)} = [time', tempX, tempY];
        lickS.len(ind(j)) = time(end) - time(1);
        lickS.conf(ind(j)) = mean(allLicks(ind(j)).conf(allLicks(ind(j)).windows(decisionID,1): allLicks(ind(j)).windows(decisionID,2)));
        
        lickDisTmpX = diff(tempX);
        lickDisTmpY = diff(tempY);
        lickSpeed = sqrt(lickDisTmpX.^2 + lickDisTmpY.^2)/(interval);
        lickVX = lickDisTmpX/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
        lickVY = lickDisTmpY/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
       
        if maxInd >= preBins+1 && maxInd < length(tempX)-2
            lickS.vXPrePeak(ind(j)) = (tempX(maxInd) - tempX(maxInd-preBins))/(preBins*interval);
            lickS.vYPrePeak(ind(j)) = (tempY(maxInd) - tempY(maxInd-preBins))/(preBins*interval);
            lickS.speedPrePeak(ind(j)) = mean(lickSpeed(maxInd-1:maxInd));
        end
        % find threshold crossing
        cross = tempY >thresh;
        cross = find(cross(1:end-1)==0 & cross(2:end)==1, 1);
        if ~isempty(cross)
            lickS.crossInd(ind(j)) = cross;
            lickS.speedAtPort(ind(j)) = lickSpeed(cross);
            lickS.vXAtPort(ind(j)) = lickVX(cross);
            lickS.vYAtPort(ind(j)) = lickVY(cross);
            lickS.crossX(ind(j)) = mean(tempX(cross));
            lickS.crossY(ind(j)) = mean(tempY(cross));
        end
        
        % use fixed inter
        if ~isempty(cross) && maxInd >= preBins+1 && maxInd < length(tempX)-2
            lickS.vXFixInter(ind(j)) = (tempX(maxInd) - tempX(cross))/((maxInd-cross)*interval);
            lickS.vYFixInter(ind(j)) = (tempY(maxInd) - tempY(cross))/((maxInd-cross)*interval);
        end

    end
end

% speedAtPort(s.allChoices==1 & ~isnan(speedAtPort)) = zscore(speedAtPort(s.allChoices==1 & ~isnan(speedAtPort)));
% speedAtPort(s.allChoices==-1 & ~isnan(speedAtPort)) = zscore(speedAtPort(s.allChoices==-1 & ~isnan(speedAtPort)));

lickS.vXFixInter(s.allChoices==1 & ~isnan(lickS.vXFixInter)) = zscore(lickS.vXFixInter(s.allChoices==1 & ~isnan(lickS.vXFixInter)));
lickS.vXFixInter(s.allChoices==-1 & ~isnan(lickS.vXFixInter)) = zscore(lickS.vXFixInter(s.allChoices==-1 & ~isnan(lickS.vXFixInter)));
lickS.vXFixInterC = lickS.vXFixInter;
lickS.vXFixInterC(s.allChoices==1) = -(lickS.vXFixInterC(s.allChoices==1));
lickS.meanXC = lickS.meanX;
lickS.meanXC(s.allChoices==1 & ~isnan(lickS.meanXC)) = zscore(lickS.meanXC(s.allChoices==1 & ~isnan(lickS.meanXC)));
lickS.meanXC(s.allChoices==-1 & ~isnan(lickS.meanXC)) = zscore(lickS.meanXC(s.allChoices==-1 & ~isnan(lickS.meanXC)));
lickS.meanXC(s.allChoices==1) = -(lickS.meanXC(s.allChoices==1));
lickS.XmaxC = lickS.XMax;
lickS.XmaxC(s.allChoices==1 & ~isnan(lickS.XmaxC)) = zscore(lickS.XmaxC(s.allChoices==1 & ~isnan(lickS.XmaxC)));
lickS.XmaxC(s.allChoices==-1 & ~isnan(lickS.XmaxC)) = zscore(lickS.XmaxC(s.allChoices==-1 & ~isnan(lickS.XmaxC)));
lickS.XmaxC(s.allChoices==1) = -(lickS.XmaxC(s.allChoices==1));
end