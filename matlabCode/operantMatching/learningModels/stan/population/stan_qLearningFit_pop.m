function [tbl, samples] = stan_qLearningFit_pop(xlFile, group, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta'});
p.addParameter('modelName', ['fourParam']);
p.addParameter('numSesh', 10);
p.addParameter('iter', 2000);
p.addParameter('warmup', []);
p.parse(varargin{:});

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

[root, sep] = currComputer();
savePath = [root 'all' sep 'allsorted' sep 'stan' sep 'population' sep p.Results.modelName sep];
if ~exist(savePath)
    mkdir(savePath);
end

switch group
    case 'dta70'
        animals = [{'CG14' 'CG15'}];
    case 'ninety'
        animals = [{'CG46' 'CG47' 'CG48' 'CG49' 'CG50' 'CG52' 'CG53' 'CG54' 'CG55' 'CG56' 'CG57' 'CG58' 'CG67'...
                    'CG68' 'CG70' 'CG75' 'CG77' 'CG78' 'CG79' 'CG80' 'CG81' 'CG82' 'CG83' 'CG84' 'CG85' 'CG86'...
                    'CG87' 'CG88' 'CG89'}];
    case 'dta90'
        animals = [{'CG84' 'CG85' 'CG86' 'CG87'}];
    case 'dtaCtrl'
        animals = [{'CG82' 'CG83' 'CG84' 'CG85'}];
end

choiceTmp = [];
outcomeTmp = [];
Tsesh = [];
for aInd = 1:length(animals)
    [~, dayList, ~] = xlsread(xlFile, animals{aInd});
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    dayList = dayList(2:end,col);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end

    for dInd = 1:p.Results.numSesh
        sessionName = dayList{dInd};
        filename = [sessionName '.asc'];
        [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
        behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

        choiceTmp = [choiceTmp {behavStruct.allChoices}];
        choiceTmp{end}(choiceTmp{end} == -1) = 0;
        outcomeTmp = [outcomeTmp {abs(behavStruct.allRewards)}]; 
        Tsesh = [Tsesh length(outcomeTmp{end})];
    end
end

M = length(animals);
N = p.Results.numSesh;
T = max(max(Tsesh));
MxN = M*N;
MxP = M*length(p.Results.paramNames);
choice = zeros(MxN, T);
outcome = zeros(MxN, T);

for m = 1:MxN
    choice(m, 1:Tsesh(m)) = choiceTmp{m};
    outcome(m, 1:Tsesh(m)) = outcomeTmp{m};
end

%create data structure to feed into stan model
session_dat = struct('M', M, 'N', N, 'T', T, 'MxN', MxN, 'MxP', MxP, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome);

%run the stan model
filePath = ['C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\population\'];
switch p.Results.modelName
    case 'fourParam'
        fit = stan('file',[filePath 'stan_qLearning_4params_pop.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'meow'
        fit = stan('file',[filePath 'stan_qLearning_4paramsMeow_pop.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_pop.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    otherwise
        error([p.Results.modelName ' does not exist'])
end

%read command line output to stall matlab until stan is finished processing
doneFlag = 0;
diary([savePath 'diaryTmp.txt']); diary off;
fid = fopen([savePath 'diaryTmp.txt'],'rt');
tmp = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
baseCount = find(~cellfun(@isempty,strfind(tmp{1}, '[100%]')) == 1);
while doneFlag == 0
    diary([savePath 'diaryTmp.txt']); pause(5); diary off;
    fid = fopen([savePath 'diaryTmp.txt'],'rt');
    tmp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    tmpCount = find(~cellfun(@isempty,strfind(tmp{1}, '[100%]')) == 1);
    if ~isempty(tmpCount)
        if length(tmpCount) == length(baseCount) + 4
            doneFlag = 1;
        end
    end
end
delete([savePath 'diaryTmp.txt'])

%extract samples from the stan fit object
samples = [];
while isempty(samples)
    samples = fit.extract('permuted',true);
end
if isfield(samples, 'y_pred')
    samples = rmfield(samples, 'y_pred');
end
%[~, tbl] = fit.print();

%generate best estimates of parameters
paramEsts = [];
for i = 1:length(p.Results.paramNames)
    paramEsts(i)  = median(eval(['samples.h_mu_' p.Results.paramNames{i}]));
end


%save the full samples
sampFile = [group category '_', p.Results.modelName];
saveFile = [sampFile '.mat'];
eval([sampFile,  ' = samples;']);
save([savePath saveFile], sampFile, 'paramEsts', 'dayList');
%save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');

%plot the distributions of the mouse-level parameters
figure; 
numParams = length(p.Results.paramNames);
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];

for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['samples.h_mu_' p.Results.paramNames{i}]) , 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    set(gca,'tickdir', 'out') 
    title(p.Results.paramNames{i})
end
titleTxt = strrep([group ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters', 'position', [-1919 1 1920 1004])
saveFigurePDF(gcf,[savePath group category '_' p.Results.modelName  '_posteriors'])
    
    
    
