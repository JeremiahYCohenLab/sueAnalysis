function compareAvgRwdRbar_opMD(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('randomSeed', 77582);
p.addParameter('modelName', 'fiveParam_rBeta_scale');
p.addParameter('revForFlag', 0);
p.parse(varargin{:});

[root,sep] = currComputer();
    
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

for i = 1:length(dayList)
    sessionName = dayList{i};
    filename = [sessionName '.asc'];
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
    seshRwdAvg(i) = sum(abs(behavStruct.allRewards)) / length(behavStruct.allRewards);


    modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
        category '_' p.Results.modelName '.mat'];
    t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, p.Results.revForFlag);
    rBarAvg(i) = mean(t.R);
    
end


figure;
scatter(seshRwdAvg, rBarAvg, 'k', 'filled')
set(gca, 'tickdir' , 'out')
xlabel('session reward rate')
ylabel('session average rBar')


