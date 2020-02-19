function [tbl, samples] = stan_spikeCountMdlFit(xlFile, sheet, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('paramNames', {'v', 'slope', 'intercept', 'sigma'});
p.addParameter('modelName', ['v']);
p.addParameter('window', [0 500]);
p.addParameter('iter', 2000);
p.addParameter('warmup', []);
p.parse(varargin{:});

tB = p.Results.window(1);
tF = p.Results.window(2);
paramInds = [1:length(p.Results.paramNames)];

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

[root, sep] = currComputer();

% load data from excel sheet
[numbers, sessionCellList, ~] = xlsread(xlFile, sheet);
revForFlag = numbers(:,1);
intanFlag = numbers(:,2);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);
if size(numbers,2) > 2
    trialList = numbers(:,3:4);
else
    trialList = nan(length(cellList),2);
end

%create savePath if doesn't already exist
savePath = [root 'spikeMdls' sep p.Results.modelName sep];
if ~exist(savePath)
    mkdir(savePath);
end


%make kernel for smoothing spikes
smoothKern = normpdf(0:5000, 0, 250);
smoothKern = smoothKern/sum(smoothKern);

sessionName = [];
numCells = length(cellList);

for currCell = 1:numCells
    if strcmp(sessionName, sessionList{currCell}) == 0
        sessionName = sessionList{currCell};
        [animalName, date] = strtok(sessionName, 'd'); 
        animalName = animalName(2:end);

        [animalName] = strtok(sessionName, 'd');
        animalName = animalName(2:end);

        if isstrprop(sessionName(end), 'alpha')
            sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
        else
            sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
        end
        sortedFolder = dir(sortedFolderLocation);


        if intanFlag(currCell)
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                [sessionData] = generateSessionData_intan_operantMatching(sessionName);
            end
        else
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                [sessionData] = generateSessionData_nL_operantMatching(sessionName);
            end
        end
        
        %get behavior data
        filename = [sessionName '.asc'];
        [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
        s = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
        
        %get all cell and session time info
        spikeFields = fields(sessionData);
        cellInds = find(~cellfun(@isempty,strfind(spikeFields, 'TT')));
        sessionTime = (sessionData(end).CSon + tF) - (sessionData(1).CSon + tB);
        csTimes = [sessionData.CSon] -  (sessionData(1).CSon + tB); 
        csTimes = [csTimes(s.responseInds) csTimes(end) + tF];  %add extra csTime for indexing purposes below
        tEndTimes = [0 csTimes(1:end-1)+2500];                  %end times of the previous trial to avoid taking spikes during previous trials
    end
    
    %get relevant behavior data
    outcomeTmp{currCell} = abs(s.allRewards); 
    Tsesh(currCell) = length(outcomeTmp{currCell});
    
    %get individual cell data and smooth
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tB);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;

    spikeTmp{currCell} = nan(1,length(s.responseInds));
    for trialInd = 1:length(s.responseInds)
        winInds = [(csTimes(trialInd) + tB + 1) : (csTimes(trialInd) + tF)];
        winInds = winInds(winInds > tEndTimes(trialInd) & winInds <  csTimes(trialInd + 1));
        if ~isempty(winInds)
            spikeTmp{currCell}(trialInd) = sum(sessionSpikes(winInds));
        end
    end
    nanInds = find(isnan(spikeTmp{currCell}));
    if ~isempty(nanInds)
        spikeTmp{currCell}(nanInds) = [];
        outcomeTmp{currCell}(nanInds) = [];
        Tsesh(currCell) = Tsesh(currCell) - length(nanInds);
    end
end

T = max(Tsesh);
N = length(cellList);
outcome = zeros(N, T);
spikes = zeros(N, T);

for i = 1:N
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
    spikes(i, 1:Tsesh(i)) = spikeTmp{i};
end

%create data structure to feed into stan model
session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'outcome', outcome, 'spikes', spikes);


%run the stan model
filePath = ['C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\functions\spikeAnalysis\spikeModels\'];
switch p.Results.modelName
    case 'v'
        fit = stan('file',[filePath 'stan_spikeCountMdl_v.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    otherwise
        fprintf([p.Results.modelName ' model does not exist'])
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
%[~, tbl] = fit.print();    
% weird error in parsing the string, which is lacking spaces

%generate best estimates of parameters for each cell
paramEsts = [];
for i = 1:length(paramInds)
    for j = 1:length(cellList)
        paramEsts(j,i)  = median(eval(['samples.' p.Results.paramNames{paramInds(i)} '(:,j)']));
    end
end


%save the full samples
sampFile = [sheet '_', p.Results.modelName];
saveFile = [sampFile '.mat'];
eval([sampFile,  ' = samples;']);
% save([savePath saveFile], sampFile, 'paramEsts', 'cellList', 'tbl');
save([savePath saveFile], sampFile, 'paramEsts', 'cellList');


%plot correlation between predicted spikes and actual
s_rho = nan(1, numCells);   s_pVal = nan(1, numCells);
p_rho = nan(1, numCells);   p_pVal = nan(1, numCells);
for i = 1:length(cellList)
    predSpikes = median(squeeze(samples.spikes_pred(:,i,:)),1);
    predSpikes = predSpikes(1:Tsesh(i));
    [s_rho(i), s_pVal(i)] = corr(predSpikes', spikes(i,1:Tsesh(i))', 'Type', 'Spearman');
    [p_rho(i), p_pVal(i)] = corr(predSpikes', spikes(i,1:Tsesh(i))', 'Type', 'Pearson');
%     mdl = fitlm(predSpikes, spikes(i,:));
%     rSqr(i) = mdl.Rsquared.Ordinary;
end


%plot the distributions of the mouse-level parameters
figure;
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numCells)', linspace(blue(2),purp(2),numCells)', linspace(blue(3),purp(3),numCells)'];

subplot(1,3,1); hold on;
for i = 1:numCells
    scatter(s_rho(i), p_rho(i), [], colors(i,:), 'filled')
end
lim = max([xlim ylim]);
xlim([0 lim]); ylim([0 lim]);
plot([0 lim], [0 lim], '--k')
xlabel('Spearmans \rho')
ylabel('Pearsons correlation coeff')
set(gca,'tickdir', 'out') 

subplot(1,3,2)
histogram(paramEsts(:,1), 15, 'Normalization', 'Probability', 'FaceColor', colors(1,:))
ylabel('Probability')
xlabel('v')
set(gca,'tickdir', 'out') 


subplot(1,3,3); hold on;
for i = 1:numCells
    scatter(paramEsts(i,1), p_rho(i), [], colors(i,:), 'filled')
end
xlabel('v')
ylabel('Pearsons correlation coeff')
set(gca,'tickdir', 'out') 

titleTxt = strrep([sheet ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')


figure;
set(gcf,'defaultAxesColorOrder',[0 0 1; 0 0 0]);
sigInds = find(p_rho > 0.25);
plotNum = ceil(length(sigInds)/2);

for i = 1:length(sigInds)
    subplot(plotNum,2,i); hold on;
    yyaxis left; plot(spikes(sigInds(i), 1:Tsesh(sigInds(i))))
    if mod(i,2) == 1
        ylabel('Spikes/s')
    end
    yL = ylim;
    yyaxis right; plotFilled([1:Tsesh(sigInds(i))], squeeze(samples.spikes_pred(:,sigInds(i),1:Tsesh(sigInds(i)))), 'k', 1)
    ylim(yL)
    if mod(i,2) == 0
        ylabel('Mdl spikes/s')
    end
    if i >= length(sigInds) - 1
        xlabel('Trials')
    end
    set(gca,'tickdir', 'out')
end