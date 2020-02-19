function [tbl, samples] = stan_spikeRateMdlFit(xlFile, sheet, varargin)

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
    choiceTmp{currCell} = s.allChoices;
    choiceTmp{currCell}(choiceTmp{currCell} == -1) = 0;
    outcomeTmp{currCell} = abs(s.allRewards); 
    Tsesh(currCell) = length(outcomeTmp{currCell});
    
    %get individual cell data and smooth
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tB);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;
    smoothSpikes = conv(sessionSpikes, smoothKern);
    smoothSpikes = smoothSpikes(1:(end-length(smoothKern)+1)) * 1000;

    spikeTmp{currCell} = nan(1,length(s.responseInds));
    for trialInd = 1:length(s.responseInds)
        winInds = [(csTimes(trialInd) + tB + 1) : (csTimes(trialInd) + tF)];
        winInds = winInds(winInds > tEndTimes(trialInd) & winInds <  csTimes(trialInd + 1));
        if ~isempty(winInds)
            spikeTmp{currCell}(trialInd) = mean(smoothSpikes(winInds));
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
choice = zeros(N, T);
spikes = zeros(N, T);

for i = 1:N
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
    choice(i, 1:Tsesh(i)) = choiceTmp{i};
    spikes(i, 1:Tsesh(i)) = spikeTmp{i};
end

%create data structure to feed into stan model
session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'outcome', outcome, 'choice', choice, 'spikes', spikes);


%run the stan model
filePath = ['C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\functions\spikeAnalysis\spikeModels\'];
switch p.Results.modelName
    case 'v'
        fit = stan('file',[filePath 'stan_spikeRateMdl_v.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'v2'
        fit = stan('file',[filePath 'stan_spikeRateMdl_v2.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'v_Rstart'
        fit = stan('file',[filePath 'stan_spikeRateMdl_v_Rstart.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sd'
        fit = stan('file',[filePath 'stan_spikeRateMdl_sd.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'vol'
        fit = stan('file',[filePath 'stan_spikeRateMdl_vol.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'peBar'
        fit = stan('file',[filePath 'stan_spikeRateMdl_peBar.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'peBar_fixedR'
        fit = stan('file',[filePath 'stan_spikeRateMdl_peBar_fixedR.stan'],'data',session_dat,'verbose',true,...
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


figure;
blue = [0 1 1];
purp = [0.7 0 1];
numParams = length(p.Results.paramNames) - 3; %plot params that are not the slope, intercept, or noise
colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];

%plot distribution of time constant parameter estimates
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    medParam = [];
    for j = 1:numCells
        medParam(j) = median(eval(['samples.' p.Results.paramNames{paramInds(i)} '(:,j)']));
    end
    histogram(medParam , 20, 'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    if i == 1
        ylabel('Probability')
    end
    set(gca,'tickdir', 'out') 
    title(p.Results.paramNames{paramInds(i)})
end
suptitle('median parameter distributions')
set(gcf,'Renderer', 'Painters') 


%plot correlation coeffs for Pearson and Spearman to check nonlinearites
figure; hold on;
colors = [linspace(blue(1),purp(1),numCells)', linspace(blue(2),purp(2),numCells)', linspace(blue(3),purp(3),numCells)'];
for i = 1:numCells
    scatter(s_rho(i), p_rho(i), [], colors(i,:), 'filled')
end
lim = max([xlim ylim]);
xlim([0 lim]); ylim([0 lim]);
plot([0 lim], [0 lim], '--k')
xlabel('Spearmans \rho')
ylabel('Pearsons correlation coeff')
set(gca,'tickdir', 'out') 

titleTxt = strrep([sheet ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')

%plot actual vs predicted firing rates by trial for significant correlations
figure;
set(gcf,'defaultAxesColorOrder',[0 0 1; 0 0 0]);
sigInds = find(s_rho > 0.25);
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
suptitle(titleTxt);

%plot scatters of actual v predicted firing rates
figure;
plotNum = ceil(length(sigInds)/4);
for i = 1:length(sigInds)
    subplot(plotNum,4,i); hold on;
    scatter(spikes(sigInds(i), 1:Tsesh(sigInds(i))), median(squeeze(samples.spikes_pred(:,sigInds(i),1:Tsesh(sigInds(i))))), ...
        [], 'k', 'filled')
    if mod(i,2) == 1
        ylabel('Mdl spikes/s')
    end
    if i >= length(sigInds) - 3
        xlabel('Spikes/s')
    end
    set(gca,'tickdir', 'out')
end
suptitle(titleTxt);


%plot predicted v residuals
figure;
plotNum = ceil(length(sigInds)/4);
for i = 1:length(sigInds)
    subplot(plotNum,4,i); hold on;
    res = spikes(sigInds(i), 1:Tsesh(sigInds(i))) - median(squeeze(samples.spikes_pred(:,sigInds(i),1:Tsesh(sigInds(i)))));
    scatter(median(squeeze(samples.spikes_pred(:,sigInds(i),1:Tsesh(sigInds(i))))), res, ...
        [], 'k', 'filled')
    if mod(i,2) == 1
        ylabel('Residuals')
    end
    if i >= length(sigInds) - 3
        xlabel('Mdl spikes/s')
    end
    set(gca,'tickdir', 'out')
end
suptitle(titleTxt);