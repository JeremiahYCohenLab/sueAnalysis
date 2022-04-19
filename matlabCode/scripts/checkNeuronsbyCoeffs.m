%% load data
load F:\tmpData\catNoInteraction.mat
catsNo = cats;
load F:\tmpData\catWithInteraction.mat
%% find neurons pos in noInter but NA in Inter
currInd = find(cats==0 & catsNo ==2);
for i = 1:length(currInd)
    spikeGLM_dF(allSessions{currInd(i)}, 'good','cellName', allUnits{currInd(i)},'regressors','1  + Qchosen + outcome*rightSide');
end
%% find neurons NA in noInter but neg in Inter
currInd = find(cats==1 & catsNo ==0);
for i = 1:length(currInd)
    spikeGLM_dF(allSessions{currInd(i)}, 'good','cellName', allUnits{currInd(i)},'regressors','1  + Qchosen + outcome*rightSide');
end  
%% check autoCorr
load F:\tmpData\catQ.mat
path = 'F:\allUnits\autoCorr\';
for i = 83:length(cats)
    if cats(i)==1
        savePath = [path 'Qneg'];
    else
        savePath = [path 'QnotNeg'];
    end
    unitAutoCorr(allSessions{i}, allUnits{i}, 'lag', 500);
    saveFigurePDF(gcf, [savePath '\' allSessions{i} '_' allUnits{i} 'autoCorr' '.pdf']);
end
%% append pdfs
savePath = ['F:\allUnits\autoCorr\QnotNeg\'];
allFiles = dir(savePath);
allFiles = {allFiles([allFiles.bytes]>0).name}';
allFiles = strcat(savePath, allFiles);
append_pdfs([savePath 'combineAutoCorrNotQneg.pdf'],allFiles{:});
    

%%