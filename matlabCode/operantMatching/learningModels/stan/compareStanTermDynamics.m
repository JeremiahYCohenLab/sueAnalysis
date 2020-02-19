function [model] = compareStanTermDynamics(xlFile, animal, pre, post, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('params', {'aN', 'pePe', 'peBar'})
p.addParameter('plotFlag', 1);
p.addParameter('bernFlag', 1);
p.addParameter('manipFlag', 1);
p.addParameter('sessionFlag', 1);
p.addParameter('numLags', 10)
p.parse(varargin{:});

[root, sep] = currComputer();

[~, dayList_pre, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList_pre, pre)) == 1);
dayList_pre = dayList_pre(2:end,col);
endInd = find(cellfun(@isempty,dayList_pre),1);
if ~isempty(endInd)
    dayList_pre = dayList_pre(1:endInd-1,:);
end

numSesh = length(dayList_pre);
numParams = length(p.Results.params);
terms_pre = cell(numSesh, numParams);
acf_pre = cell(numParams);
acf_post = cell(numParams);

meow = [];

for currS = 1:numSesh
    sessionName = dayList_pre{currS};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    
    if p.Results.manipFlag
        modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep ...
            p.Results.modelName '_manip' sep animalName pre post '_' p.Results.modelName '_manip.mat'];
        t = generateStanModelTermsManip(p.Results.modelName, modelPath, sessionName, 1, p.Results.sessionFlag);
    else
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            pre '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            pre '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, p.Results.sessionFlag);
    end
    
    for currP = 1:numParams
        tmp = eval(['t.' p.Results.params{currP}]);
        terms_pre{currS, currP} = tmp;
        acf_pre{currP} = ([acf_pre{currP}; autocorr(tmp, p.Results.numLags)']);
    end

    meow = [meow; t.aN t.o];

end

termsMatx_pre = [];
for currS = 1:numSesh
    termsMatx_pre = [termsMatx_pre terms_pre{currS, 1}'];
end
if numParams > 1
    trials = length(termsMatx_pre);
    for currP = 2:numParams
        tmp = [];
        for currS = 1:numSesh
            tmp = [tmp terms_pre{currS, currP}'];
        end
        termsMatx_pre = [termsMatx_pre; tmp];
    end
end



[~, dayList_post, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList_post, post)) == 1);
dayList_post = dayList_post(2:end,col);
endInd = find(cellfun(@isempty,dayList_post),1);
if ~isempty(endInd)
    dayList_post = dayList_post(1:endInd-1,:);
end


numSesh = length(dayList_post);
terms_post = cell(numSesh, numParams);

for currS = 1:numSesh
    sessionName = dayList_post{currS};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    
    if p.Results.manipFlag
        modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep ...
            p.Results.modelName '_manip' sep animalName pre post '_' p.Results.modelName '_manip.mat'];
        t = generateStanModelTermsManip(p.Results.modelName, modelPath, sessionName, 0, p.Results.sessionFlag);
    else
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            post '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            post '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, p.Results.sessionFlag);
    end
    
    for currP = 1:numParams
        tmp = eval(['t.' p.Results.params{currP}]);
        terms_post{currS, currP} = tmp;
        acf_post{currP} = ([acf_post{currP}; autocorr(tmp, p.Results.numLags)']);
    end
    
    meow = [meow; t.aN t.o];
    
end

termsMatx_post = [];
for currS = 1:numSesh
    termsMatx_post = [termsMatx_post terms_post{currS, 1}'];
end
if numParams > 1
    trials = length(termsMatx_post);
    for currP = 2:numParams
        tmp = [];
        for currS = 1:numSesh
            tmp = [tmp terms_post{currS, currP}'];
        end
        termsMatx_post = [termsMatx_post; tmp];
    end
end

colors = cool(numParams*2);
figure;
for currP = 1:numParams
    subplot(2,numParams,currP); hold on;
    histogram(termsMatx_pre(currP,:), 'Normalization', 'Probability', 'FaceColor', colors(currP,:))
    histogram(termsMatx_post(currP,:), 'Normalization', 'Probability', 'FaceColor', colors(currP+numParams,:))
    xlabel('term value')
    ylabel('probability')
    legend('pre', 'post')
    set(gca, 'tickdir', 'out')
    title(p.Results.params{currP})
    
    subplot(2, numParams, numParams+currP); hold on;
    plotFilled([0:p.Results.numLags], acf_pre{currP}, colors(currP,:));
    plotFilled([0:p.Results.numLags], acf_post{currP}, colors(currP+numParams,:));
    xlabel('lag')
    ylabel('correlation coeff')
    set(gca, 'tickdir', 'out')
end

titleTxt = strrep([animal ' ' p.Results.modelName], '_', ' ');
suptitle(titleTxt)
set(gcf, 'renderer', 'painters', 'position', [-1516 284 890 667])








