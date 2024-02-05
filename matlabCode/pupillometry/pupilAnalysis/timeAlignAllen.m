function [errorProp,csFT, qualInd, ratioMax] = timeAlignAllen(session, plotFlag, saveFlag)
%% default time window is -2s to 10s to cue time
% calculate time projection of behavior time onto pupil time 
unblockedError = false; 
iterMaxReach = false;
errorRate = 0.05;
errorThresh = 6; %no. of frames allowed for mis-alignment
iter = 0;
errorProp = NaN;
csFT = NaN;
ratioMax = NaN;
ledLLThresh = 0.99999;
positionThresh = 0.999999999999;
ledDisThresh = 300;
%% training set
[root,sep] = currComputer;
   
%% behavior
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];

if isstrprop(session(end), 'alpha')
    behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_sessionData_behav.mat'];
else
    behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_behav.mat'];
end


if exist(behSessionDataPath,'file')
    load(behSessionDataPath)
else
%     behSessionDate = generateSessionData_behav_operantMatching_RwdDelay(session);
    behSessionData = generateSessionData_operantMatchingDecoupledRwdDelay(session);
end


if isstrprop(session(end), 'alpha')
    savepath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep];
else
    savepath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep];
end

iterMax = 0.05*length(behSessionData);
%% load diameter and position
videopath = [root animalName sep sessionFolder sep 'pupil'];
list = dir(videopath);
expression = ['^' session 'DLC' '\w*' '200000.csv' '$'];
expressionSkeleton = ['^' session 'DLC' '\w*' 'shuffle1_200000_skeleton.csv' '$'];
if isempty(find(~cellfun(@isempty, cellfun(@(x) regexp(x, expression), {list.name}, 'UniformOutput', false)), 1))
   fprintf([session ' no pupil video \n'])
   return
end
position = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expression), {list.name}, 'UniformOutput', false))).name;
if isempty(position)
    fprintf([session ' pupil video not amalyzed \n'])
    return
end
positionRaw = csvread([videopath sep position],3,0);

%% Align time by finding maximum projection
ratio =[22:0.01:23];
%ledLL = sign(positionRaw(:,16).*positionRaw(:,10) - ledLLThresh);
ledLL = sign(positionRaw(:,16) - ledLLThresh);
%% led position
ledx = mean(positionRaw(positionRaw(:,16)>positionThresh,14));
ledy = mean(positionRaw(positionRaw(:,16)>positionThresh,15));
dis = sqrt((positionRaw(:,14)-ledx).^2 + (positionRaw(:,15)-ledy).^2);
ledLL(dis > ledDisThresh) = -1; % exlude led too far away. (likely to be false positive)
ledQual = positionRaw(:,10);
%ledLL = ledLL .* (ledQual.^2);
csT = [behSessionData(cellfun(@(x) strcmp(x,'CSplus'), {behSessionData.trialType})).CSon];
csT = csT - csT(1);
filter = ones(1,round(0.5 * mean(ratio)));
accLL = conv(filter, 0.5 * (ledLL + 1));
accLL = accLL(length(filter) : end);
startFrameTemp = find(accLL > 0.7 * length(filter) & ledLL > 0 );
startFrameTemp = startFrameTemp(ledLL(startFrameTemp-1)<0);
startFrame = startFrameTemp(1);
startTrial = 1;
[ratioMax, ~, ~] = timeProjectionOpti2(ledLL, csT(1:50), ratio, startFrame, startTrial);
%% check first error
cskernel = ones(1, round(0.5 * ratioMax));
[~, csF, csFT] = timeProjection(ledLL, csT, startFrame, startTrial, ratioMax);
cs = zeros(1,length(csT));
csQu = zeros(1,length(csT));
for j = 1:length(csT)
    if j < length(csFT)
        cs(j) = csF(csFT(j):csFT(j) + round(0.5 * ratioMax))*ledLL(csFT(j):csFT(j) + round(0.5 * ratioMax));
        csQu(j) = mean(ledQual(csFT(j):csFT(j) + round(0.5 * ratioMax)));
    else
        cs(j) = 0;
        csQu(j) = 0;
    end
    
end
%%
cs = length(cskernel) - csQu .* (length(cskernel) - cs);  
error = find(cs < length(cskernel) - errorThresh);
%% updating by adding start points
while length(error) > errorRate*length(csT) && iter < iterMax
    errorGaps = diff(error);% gaps between error trials
    errorSeq = find(errorGaps(1:end-1) == 1 & errorGaps(2:end)== 1);
    % break if errors are not blocked for certain number, indicating
    % shifting, here is 3.
    if isempty(errorSeq)
        unblockedError = true; 
        break % break if errors do not form cluster
    end
    % separating point
    for a = 1:length(errorSeq)
        if (sum(ledLL(csFT(error(errorSeq(a))) : csFT(error(errorSeq(a))) + 4)) > 0) && min(abs(error(errorSeq(a))+1 - startTrial))>5
            errorSeq = error(errorSeq(a));
            break;
        else
            if a == length(errorSeq)
              errorSeq = error(errorSeq(1));
            end
        end
    end
    
    % adding starting points
    startFramePlus = startFrameTemp(startFrameTemp > csFT(errorSeq) + 50);
    if startFramePlus
        startFramePlus = startFramePlus(1);
    else
        break
    end
    startFrame = sort([startFrame startFramePlus]);
    startTrial = sort([startTrial errorSeq + 1]);
    % optimizing with new startpoints
    [ratioMax, ~, ~] = timeProjectionOpti2(ledLL, csT, ratio, startFrame, startTrial);   

    % calculate error number
    cskernel = ones(1, round(0.5* ratioMax));
    [~, csF, csFT] = timeProjection(ledLL, csT, startFrame, startTrial, ratioMax);
    
    for j = 1:length(csT)
        if j < length(csFT)
            cs(j) = csF(csFT(j):csFT(j) + round(0.5 * ratioMax))*ledLL(csFT(j):csFT(j) + round(0.5 * ratioMax));
            csQu(j) = mean(ledQual(csFT(j):csFT(j) + round(0.5 * ratioMax)));
        else
            cs(j) = 0;
            csQu(j) = 0;
        end
    end
    cs = length(cskernel) - csQu .* (length(cskernel) - cs); 
    error = find(cs < length(cskernel) - errorThresh);
    iter = iter + 1;
end
pupilIdx = ~(cs < length(cskernel) - errorThresh); %% well aligned ones
if length(error) > errorRate*length(csT) && iter == iterMax
    iterMaxReach = true;
end
pupilInds = pupilIdx==1;
errorProp = length(error)/length(csT);
%% plotting everything
[~, p] = timeProjectionOpti(ledLL, csT, startFrame, startTrial, ratio); 
if plotFlag
figure; hold on;
screenSize = get(0,'Screensize');
screenSize(4) = screenSize(4) - 100;
set(gcf, 'Position', screenSize)
sgtitle(session)
subplot(1,3,1)
plot(ratio, p);
line([ratioMax ratioMax], minmax(p), 'Color', 'red');
text(mean(ratio), max(p)+0.3, sprintf('errorRate %.2g', errorProp));
if unblockedError
    text(mean(ratio), max(p)+0.5, 'unblockedErrors');
end
if iterMaxReach
    text(mean(ratio), max(p)+0.5, 'maxIter reached');
end
xlim([min(ratio) max(ratio)])

subplot(2,3,[2,3]); hold on;
plot(csF, 'b');
plot(ledLL, 'r');
for j=1:length(csFT)
    text(csFT(j), 1, num2str(j), 'FontSize', 10);
end
scatter(csFT(error) + 0.5*(length(cskernel)), 0.5*ones(size(error)), 30, 'filled')
scatter(startFrame, 0.75*ones(size(startFrame)), 30, 'red')
xlim([1 2000]);
ylim([0 1.2]);

subplot(2,3,[5,6]); hold on;
plot(csF, 'b');
plot(ledLL, 'r');
for j=1:length(csFT)
    text(csFT(j), 1, num2str(j), 'FontSize', 10);
end
scatter(csFT(error) + 0.5*(length(cskernel)), 0.5*ones(size(error)), 30, 'filled')
scatter(startFrame, 0.75*ones(size(startFrame)), 30, 'red')
xlim([csFT(end)-2000 csFT(end)]);
ylim([0 1.2]);
end

%% calculate frames correspond to all cues on
csminusInds = find(cellfun(@(x) strcmp(x,'CSminus'), {behSessionData.trialType})>0);
csplusInds = find(cellfun(@(x) strcmp(x,'CSplus'), {behSessionData.trialType})>0);
% csplus
cueFT = zeros(1,length(behSessionData));
cueFT(csplusInds(1:length(csFT)))= csFT;
% csminus
csminusInfer = zeros(1,length(csminusInds));
for k = 1:length(csminusInds)
    p = 1;
    if csminusInds(k)==1
        csminusInfer = zeros(1,length(csminusInds)-1);
        continue
    end       
    while strcmp(behSessionData(csminusInds(k)-p).trialType, 'CSminus')
        p = p+1;
    end
     csminusInfer(k)=csminusInds(k)-p;
    cueFT(csminusInds(k)) = cueFT(csminusInfer(k)) + round(ratioMax*0.001*(behSessionData(csminusInds(k)).CSon - behSessionData(csminusInds(k)-p).CSon));
end
%% calculate frame quality for all cues
% create good index in csInds by pupilInds, diameterInds
qualInd = zeros(1,length(behSessionData));
qualF = positionRaw(:,4) > 0.95 &  positionRaw(:,7) > 0.95; 
for j = 1:length(qualInd)
    if cueFT(j)
        minF = max(1, min(cueFT(j) - 30, length(qualF)));
        maxF = min(cueFT(j) + 100, length(qualF));
        qualInd(j) = sum(qualF(minF : maxF));% frame quality around the time of cue
    end
end
qualInd = find(qualInd > 100); % 100 out of 130 frames are good
pupilInds = csplusInds(pupilInds);% convert from csplus indexing to all trial type indexing 
qualplus = intersect(qualInd, pupilInds);% both good frame quality and good alignment
lowqualCSplus = csplusInds(~pupilIdx);
poorMinus = ismember(csminusInfer,lowqualCSplus);
qualminus = intersect(qualInd, csminusInds(~poorMinus));% good frame quality and good alighment of previos csplus
qualInd = sort([qualplus, qualminus]);
qualInd = ismember(1:length(cueFT),qualInd);
FR = ratioMax;

%% load diameter
% load diameter
skeleton = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expressionSkeleton), {list.name}, 'UniformOutput', false))).name;
diaRaw = csvread([videopath sep skeleton], 2, 0);
qualF(end) = 1; % make sure not having a lot of NaN
dia = interp1(find(qualF>0),diaRaw(qualF,2),1:length(diaRaw(:,2)));
% filter
fc = 5;
fs = FR;
[b,a] = butter(2,fc/(fs/2),'low');
dia = filtfilt(b,a,dia);
% realign to continuous time
% put realigned pupil frame back to linear time
diaRealign = NaN(1,cueFT(1)-1 + ceil(FR/1000 * (behSessionData(end).CSon + 10000 - behSessionData(1).CSon)));
diaRealign(1:cueFT(1)-1) = dia(1:cueFT(1)-1); % quality before cue is usually good
for j = 1:length(cueFT)
    if qualInd(j)
        startF = cueFT(1) + round(FR/1000*(behSessionData(j).CSon - behSessionData(1).CSon));
        if j~=length(cueFT)
            endF = cueFT(1) + round(FR/1000*(behSessionData(j+1).CSon - behSessionData(1).CSon));
        else
            endF = min(cueFT(1) + round(FR/1000*(behSessionData(j).CSon + 10000 - behSessionData(1).CSon)), length(diaRealign));
        end

        if cueFT(j)+endF-startF > length(dia) % in case pupil ended early
            endF = length(dia)+startF-cueFT(j);
        end
        diaRealign(startF:endF) = dia(cueFT(j):cueFT(j)+endF-startF); 
    end      
end

% zscoring
diaZ = dia;
diaZ(~isnan(dia)) = zscore(dia(~isnan(dia)));

% put aligned pupil diameter together
preLen = round(2*FR);
postLen = round(10*FR);
sessionPupilCue = nan(length(cueFT),preLen+postLen+1);
sessionPupilCueZ = nan(length(cueFT),preLen+postLen+1);
for i = 1:length(cueFT)
    if cueFT(i)~=0 && cueFT(i)<=length(dia)
        startF = max([1, cueFT(i)-preLen]);
        if i==length(cueFT)
            endF = min([length(ledLL), cueFT(i)+postLen]);
        else
            endF = min([cueFT(i)+postLen, cueFT(i+1)]);
        end

        sessionPupilCue(i, preLen+1-(cueFT(i)-startF):preLen+1) = dia(startF:cueFT(i));
        sessionPupilCue(i, preLen+2:preLen+1+endF-cueFT(i)) = dia(cueFT(i)+1:endF);

        sessionPupilCueZ(i, preLen+1-(cueFT(i)-startF):preLen+1) = diaZ(startF:cueFT(i));
        sessionPupilCueZ(i, preLen+2:preLen+1+endF-cueFT(i)) = diaZ(cueFT(i)+1:endF);
    end
end

responseInds = find(~isnan([behSessionData.rewardTime]));
lickLat = [behSessionData(responseInds).respondTime] - [behSessionData(responseInds).CSon];
choiceFT = cueFT(responseInds) + round(FR*lickLat/1000);

sessionPupilChoice = nan(length(choiceFT),preLen+postLen+1);
sessionPupilChoiceZ = nan(length(choiceFT),preLen+postLen+1);
for i = 1:length(choiceFT)
    if cueFT(i)> 0 && cueFT(i)<=length(dia)
        startF = max([1, choiceFT(i)-preLen]);
        if i==length(choiceFT)
            endF = min([length(ledLL), choiceFT(i)+postLen]);
        else
            endF = min([choiceFT(i)+postLen, choiceFT(i+1)]);
        end

        sessionPupilChoice(i, preLen+1-(choiceFT(i)-startF):preLen+1) = dia(startF:choiceFT(i));
        sessionPupilChoice(i, preLen+2:preLen+1+endF-choiceFT(i)) = dia(choiceFT(i)+1:endF);

        sessionPupilChoiceZ(i, preLen+1-(choiceFT(i)-startF):preLen+1) = diaZ(startF:choiceFT(i));
        sessionPupilChoiceZ(i, preLen+2:preLen+1+endF-choiceFT(i)) = diaZ(choiceFT(i)+1:endF);
    end
end

if saveFlag
    if isempty(dir(savepath))
        mkdir(savepath)
    end
    eyeBlurryLate = 0;
    save([savepath session '_pupil.mat'], 'eyeBlurryLate', 'sessionPupilCue','sessionPupilCueZ', 'sessionPupilChoice','sessionPupilChoiceZ','qualInd', 'cueFT', 'FR','iter','errorProp', 'dia', 'diaRealign');
end





