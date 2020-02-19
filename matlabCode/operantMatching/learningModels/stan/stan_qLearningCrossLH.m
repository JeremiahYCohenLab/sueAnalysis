function [LH] = stan_qLearningCrossLH(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('fixedParams', 0);
p.addParameter('params', ['fiveParamOfourStart_preS']);
p.addParameter('paramInds', [1 2 3 4 5]);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'v'});
p.addParameter('animalName', []);
p.addParameter('modelType', ['fourParam']);
p.parse(varargin{:});

[root, sep] = currComputer();

[~, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

for i = 1:length(dayList)
    sessionName = dayList{i};
    filename = [sessionName '.asc'];
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

    choiceTmp{i} = behavStruct.allChoices;
    choiceTmp{i}(choiceTmp{i} == -1) = 2;
    outcomeTmp{i} = abs(behavStruct.allRewards); 
    Tsesh(i,1) = length(outcomeTmp{i});
end


inds = [1:length(dayList)];
aInds = inds(1:2:end);      bInds = inds(2:2:end);
TseshA = Tsesh(aInds);      TseshB = Tsesh(bInds);
Ta = max(TseshA);           Tb = max(TseshB);
Na = length(aInds);         Nb = length(bInds);
choiceA = zeros(Na, Ta);     choiceB = zeros(Nb, Tb);
outcomeA = zeros(Na, Ta);    outcomeB = zeros(Nb, Tb);

for i = 1:length(aInds)
    choiceA(i, 1:Tsesh(aInds(i))) = choiceTmp{aInds(i)};
    outcomeA(i, 1:Tsesh(aInds(i))) = outcomeTmp{aInds(i)};
end
for i = 1:length(bInds)
    choiceB(i, 1:Tsesh(bInds(i))) = choiceTmp{bInds(i)};
    outcomeB(i, 1:Tsesh(bInds(i))) = outcomeTmp{bInds(i)};
end

session_dat = struct('N',Na,'T',Ta, 'Tsesh', TseshA, 'choice', choiceA, 'outcome', outcomeA);

filePath = ['C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\'];
switch p.Results.modelType
    case 'fourParam'
        fit = stan('file',[filePath 'stan_qLearning_4params.stan'],'data',session_dat,'verbose',true);
    case 'fiveParamO'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency.stan'],'data',session_dat,'verbose', true);
    case 'fiveParamO_rBarStart'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_rBarStart.stan'],'data',session_dat,'verbose', true);
    case 'fiveParamO_strongPriors'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_strongPriors.stan'],'data',session_dat,'verbose', true);
    case 'fiveParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta.stan'],'data',session_dat,'verbose', true);
    case 'sixParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_peBeta.stan'],'data',session_dat,'verbose', true);
    otherwise
        fprintf([p.Results.modelType ' model does not exist'])
end


doneFlag = 0;
doneCount = 0;
while doneFlag == 0
    diary('diaryTmp.txt'); pause(30); diary off;
    fid = fopen('diaryTmp.txt','rt');
    tmp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    tmpCount = find(~cellfun(@isempty,strfind(tmp{1}, '[100%]')) == 1);
    if ~isempty(tmp)
        doneCount = doneCount + length(tmpCount)
    end
    if doneCount == 4
        doneFlag = 1;
    end
end

samplesA = fit.extract('permuted',true);

paramEstsA = [];
for i = 1:length(p.Results.paramInds)
    [binNum, edges]  = discretize(eval(['samplesA.mu_' p.Results.paramNames{p.Results.paramInds(i)}]),100);
    binInd = mode(binNum);
    paramEstsA(i) = mean([edges(binInd) edges(binInd + 1)]);
end

for i = 1:length(aInds)
    [mdl] = calculateChoiceProb_opMD(dayList{aInds(i)}, 'modelNames', {p.Results.modelType}, 'params', {paramEstsA})
    LHa(i) = eval(['mdl.' p.Results.modelType '.LH']);
end


switch p.Results.modelType
    case 'fourParam'
        fit = stan('file',[filePath 'stan_qLearning_4params.stan'],'data',session_dat,'verbose',true);
    case 'fiveParamO'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency.stan'],'data',session_dat,'verbose', true);
    case 'fiveParamO_rBarStart'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_rBarStart.stan'],'data',session_dat,'verbose', true);
    case 'fiveParamO_strongPriors'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_strongPriors.stan'],'data',session_dat,'verbose', true);
    case 'fiveParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta.stan'],'data',session_dat,'verbose', true);
    case 'sixParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_peBeta.stan'],'data',session_dat,'verbose', true);
    otherwise
        fprintf([p.Results.modelType ' model does not exist'])
end

doneFlag = 0;
doneCount = 0;
while doneFlag == 0
    diary('diaryTmp.txt'); pause(30); diary off;
    fid = fopen('diaryTmp.txt','rt');
    tmp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    tmpCount = find(~cellfun(@isempty,strfind(tmp{1}, '[100%]')) == 1);
    if ~isempty(tmp)
        doneCount = doneCount + length(tmpCount)
    end
    if doneCount == 4
        doneFlag = 1;
    end
end

samplesB = fit.extract('permuted',true);

paramEstsB = [];
for i = 1:length(p.Results.paramInds)
    [binNum, edges]  = discretize(eval(['samplesA.mu_' p.Results.paramNames{p.Results.paramInds(i)}]),100);
    binInd = mode(binNum);
    paramEstsB(i) = mean([edges(binInd) edges(binInd + 1)]);
end

for i = 1:length(bInds)
    [mdl] = calculateChoiceProb_opMD(dayList{bInds(i)}, 'modelNames', p.Results.modelType, 'params', paramEstsA)
    LHb(i) = eval(['mdl.' p.Results.modelType 'LH']);
end

LH = sum([LHa LHb]) / sum(Tsesh);


end



