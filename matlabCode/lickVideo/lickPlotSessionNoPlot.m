function kine = lickPlotSessionNoPlot(session)
    [root, sep] = currComputer();
    pd = parseSessionString_df(session, root, sep);
    savepath = [pd.sortedFolder session '_tongue.mat'];
    if exist(savepath, "file")
        load(savepath);
    else
        fprintf([session ' no tongue. \n'])
        return
    end
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
%%
ind  = 1:length(allLicks);
meanX = NaN(1, length(ind));
meanY = NaN(1, length(ind));
meanXPre = NaN(1, length(ind));
meanXMax = NaN(1, length(ind));
meanXMaxZ = NaN(1, length(ind));
meanYMax = NaN(1, length(ind));
meanYMaxZ = NaN(1, length(ind));
maxS = NaN(1, length(ind));
speedAtPort = NaN(1, length(ind));
vXAtPort = NaN(1, length(ind));
vYAtPort = NaN(1, length(ind));
vXPrePeak = NaN(1, length(ind));
vYPrePeak = NaN(1, length(ind));
vXFixInter = NaN(1, length(ind));
vYFixInter = NaN(1, length(ind));
speedFixInter = NaN(1, length(ind));
speedPrePeak = NaN(1, length(ind));
crossInd = NaN(1, length(ind));
maxIndAll = NaN(1, length(ind));
crossX = NaN(1, length(ind));
crossY = NaN(1, length(ind));
allX = cell(length(ind), 1);
allY = cell(length(ind), 1);
thresh = 0;
preBins = 6;

for j = 1:length(ind)
    if allLicks(ind(j)).decisionID
        decisionID = allLicks(ind(j)).decisionID;
        interval = allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1);
        disRL = norm(allLicks(ind(j)).portL - allLicks(ind(j)).portR);
        temp = allLicks(ind(j)).allLicksFilt{decisionID};
        tempX = (temp(:,1) - allLicks(ind(j)).midPoint(1)) * 4/disRL;
        tempY = -(temp(:,2) - allLicks(ind(j)).midPoint(2)) * 4/disRL;
        allX{j} = tempX;
        allY{j} = tempY;
        [meanYMax(ind(j)), maxInd] = max(tempY);
        maxIndAll(ind(j)) = maxInd;
        meanXMax(ind(j)) = tempX(maxInd);
        meanXPre(ind(j)) = mean(tempX(1:maxInd));
        
        meanX(ind(j)) = mean(tempX);
        meanY(ind(j)) = mean(tempY);
        time = allLicks(ind(j)).time(allLicks(ind(j)).windows(decisionID,1): allLicks(ind(j)).windows(decisionID,2));

        lickDisTmpX = diff(tempX);
        lickDisTmpY = diff(tempY);
        lickSpeed = sqrt(lickDisTmpX.^2 + lickDisTmpY.^2)/(interval);
        maxS(ind(j)) = max(lickSpeed);
        lickVX = lickDisTmpX/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
        lickVY = lickDisTmpY/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
       
        if maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXPrePeak(ind(j)) = (tempX(maxInd) - tempX(maxInd-preBins))/(preBins*interval);
            vYPrePeak(ind(j)) = (tempY(maxInd) - tempY(maxInd-preBins))/(preBins*interval);
            speedPrePeak(ind(j)) = mean(lickSpeed(maxInd-1:maxInd));
        end
        % find threshold crossing
        cross = tempY >thresh;
        cross = find(cross(1:end-1)==0 & cross(2:end)==1, 1);
        if ~isempty(cross)
            crossInd(ind(j)) = cross;
            speedAtPort(ind(j)) = lickSpeed(cross);
            vXAtPort(ind(j)) = lickVX(cross);
            vYAtPort(ind(j)) = lickVY(cross);
            crossX(ind(j)) = mean(tempX(cross));
            crossY(ind(j)) = mean(tempY(cross));
        end
        
        % use fixed inter
        if ~isempty(cross) && maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXFixInter(ind(j)) = (tempX(maxInd) - tempX(cross))/((maxInd-cross)*interval);
            vYFixInter(ind(j)) = (tempY(maxInd) - tempY(cross))/((maxInd-cross)*interval);
            speedFixInter(ind(j)) = norm([vXFixInter(ind(j)), vYFixInter(ind(j))]);
        end
    else
        % fprintf([num2str(j) '\n']);
    end
end
%% regularization
vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)));
vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)));
vXFixInterAbs = vXFixInter;
vXFixInterAbs(s.allChoices==1) = -vXFixInterAbs(s.allChoices==1);

meanXMax(~isnan(meanXMax) & s.allChoices==1) = meanXMax(~isnan(meanXMax) & s.allChoices==1)/abs(mean(meanXMax(~isnan(meanXMax) & s.allChoices==1), 'all', 'omitnan'));
meanXMax(~isnan(meanXMax) & s.allChoices==-1) = meanXMax(~isnan(meanXMax) & s.allChoices==-1)/abs(mean(meanXMax(~isnan(meanXMax) & s.allChoices==-1), 'all', 'omitnan'));


meanXMaxZ(~isnan(meanXMax) & s.allChoices==1) = zscore(meanXMax(~isnan(meanXMax) & s.allChoices==1));
meanXMaxZ(~isnan(meanXMax) & s.allChoices==-1) = zscore(meanXMax(~isnan(meanXMax) & s.allChoices==-1));

meanXMaxC = meanXMaxZ;
meanXMaxC(s.allChoices==1) = - meanXMaxC(s.allChoices==1);

meanYMax(~isnan(meanYMax) & s.allChoices==1) = meanYMax(~isnan(meanYMax) & s.allChoices==1)/abs(mean(meanYMax(~isnan(meanYMax) & s.allChoices==1), 'all', 'omitnan'));
meanYMax(~isnan(meanYMax) & s.allChoices==-1) = meanYMax(~isnan(meanYMax) & s.allChoices==-1)/abs(mean(meanYMax(~isnan(meanYMax) & s.allChoices==-1), 'all', 'omitnan'));

meanYMaxZ(~isnan(meanYMax) & s.allChoices==1) = zscore(meanYMax(~isnan(meanYMax) & s.allChoices==1));
meanYMaxZ(~isnan(meanYMax) & s.allChoices==-1) = zscore(meanYMax(~isnan(meanYMax) & s.allChoices==-1));

vYFixInter(s.allChoices==1 & ~isnan(vYFixInter)) = zscore(vYFixInter(s.allChoices==1 & ~isnan(vYFixInter)));
vYFixInter(s.allChoices==-1 & ~isnan(vYFixInter)) = zscore(vYFixInter(s.allChoices==-1 & ~isnan(vYFixInter)));

speedFixInter(s.allChoices==1 & ~isnan(speedFixInter)) = zscore(speedFixInter(s.allChoices==1 & ~isnan(speedFixInter)));
speedFixInter(s.allChoices==-1 & ~isnan(speedFixInter)) = zscore(speedFixInter(s.allChoices==-1 & ~isnan(speedFixInter)));

maxS(s.allChoices==1 & ~isnan(maxS)) = zscore(maxS(s.allChoices==1 & ~isnan(maxS)));
maxS(s.allChoices==-1 & ~isnan(maxS)) = zscore(maxS(s.allChoices==-1 & ~isnan(maxS)));

meanXMax = -meanXMax;




%% combine all data
kine.xMax = meanXMax;
kine.xMaxC = meanXMaxC;
kine.yMax = meanYMax;
kine.yMaxZ = meanYMaxZ;
kine.vXFixInter = vXFixInter;
kine.vYFixInter = vYFixInter;
kine.speedFixInter = speedFixInter;
kine.vXFixInterC = vXFixInterAbs;
kine.allX = allX;
kine.allY = allY;
kine.speedPrePeak = speedPrePeak;
kine.meanX = meanX;
kine.meanY = meanY;
kine.maxS = maxS;
kine.maxInd = maxIndAll;

end 