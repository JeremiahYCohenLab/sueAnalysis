function [latencyToResponse] = optoIDrawTrace(filename, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('filtFlag', 0)
p.addParameter('delayComp',0)
p.addParameter('TTdelay', 0)
p.addParameter('CSCdelay', 0)
p.addParameter('spikeDelay', 0)
p.addParameter('intanFlag', 0)
p.addParameter('cellName', [])
p.addParameter('plotFlag', 1)
p.parse(varargin{:});


[root, sep] = currComputer();

[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'opto'];

if p.Results.filtFlag == 1
    filtDataPath = [root animalName sep sessionFolder sep 'ephys' sep 'opto' sep 'filtered' sep];
end

unsortedDataPath = [root animalName sep sessionFolder sep 'ephys' sep 'opto' sep];
if ~exist(unsortedDataPath)
    latencyToResponse = 0;
    return;
end

if exist(strcat(sessionDataPath, sep, 'laserStruct.mat'),'file')
    load(strcat(sessionDataPath, sep, 'laserStruct.mat'))
else
    [laserStruct] = findLaserEvents(filename);
    load(strcat(sessionDataPath, sep, 'laserStruct.mat'))
end


if p.Results.delayComp == 1
    p.Results.spikeDelay = p.Results.TTdelay - p.Results.CSCdelay;
end

% find cluster data
spikeFields = fields(laserStruct);
if isempty(p.Results.cellName)
    clust = find(~cellfun(@isempty,strfind(spikeFields,'SS')) & ~cellfun(@isempty,strfind(spikeFields,'TT')));
else
    clust = find(~cellfun(@isempty,strfind(spikeFields,p.Results.cellName)));
    if isempty(clust)
        latencyToResponse = 0;
        return;
    end
end


sCon = 0.000000030518510385491027*1e6; %converts spike amplitude to uV
tB = 500; % ms before first pulse
tF = 500; % ms after last pulse
numPulses = 10; %number of pulses in each train
numTrains = length(laserStruct.laserOn)/numPulses;
rasterLength = length(-1*tB:laserStruct.laserOff(numPulses)-laserStruct.laserOn(1)+tF);

saveDir = [root animalName sep sessionFolder sep 'figures'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

for i = clust'
%     CSCchann = spikeFields{i}(strfind(spikeFields{i},'CSC'):end);
    
    TTchann = str2double(spikeFields{i}(3));
    CSCchann{1} = ['CSC' num2str(4*TTchann-3)];
    CSCchann{2} = ['CSC' num2str(4*TTchann-2)]; 
    CSCchann{3} = ['CSC' num2str(4*TTchann-1)];
    CSCchann{4} = ['CSC' num2str(4*TTchann-0)];
    trainInd = 1:numPulses:length(laserStruct.laserOn);
    
    latencyToResponse = []; % latency of light-evoked responses
    numResponses = []; % success of light-evoked response
    spikeTimes_light_TT = []; % light-evoked spike in TT time
    spikeTimes_light_CSC = []; % light-evoked spike in CSC time; always the first index of the 1x512 spike vector
    spikeCSC_light = []; % light-evoked spike waveform; 1x512 spike vector
    spikeInd_light = []; % index in the 1x512 spike vector that corresponds to spike

    spikeTimes_spont_TT = [];
    spikeTimes_spont_CSC = [];
    spikeCSC_spont = [];
    spikeInd_spont = [];
    
    for j = trainInd
        lL = laserStruct.laserOn(j) - tB;
        uL = laserStruct.laserOff(j+numPulses-1) + tF;
        spikeInds = laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i}) >= lL & laserStruct.(spikeFields{i}) < uL) - laserStruct.laserOn(j) + tB;
        spikeRast(ceil(j/numPulses),1) = {spikeInds};
        
        laserOnInds = laserStruct.laserOn(laserStruct.laserOn >= lL & laserStruct.laserOn < uL) - laserStruct.laserOn(j) + tB;
        laserOffInds = laserStruct.laserOff(laserStruct.laserOff >= lL & laserStruct.laserOff < uL) - laserStruct.laserOn(j) + tB;
        laserRast(ceil(j/numPulses),1) = {sort([laserOnInds laserOffInds])};
        
        for k = 1:numPulses
            tempSpike = laserStruct.(spikeFields{i})(find(laserStruct.(spikeFields{i}) > laserStruct.laserOn(j+k-1),1));
            tempLat = tempSpike - laserStruct.laserOn(j+k-1);
            if tempLat < 30 %30ms limit to evoke a spike
                latencyToResponse(ceil(j/numPulses),k) = tempLat - p.Results.TTdelay/1000;
                numResponses(ceil(j/numPulses),k) = 1;
%                 spikeTimes_light_TT(ceil(j/numPulses),k) = tempSpike; %(tempSpike*1000)-spikeDelay;
                spikeTimes_light_TT = [spikeTimes_light_TT tempSpike];
            else
                latencyToResponse(ceil(j/numPulses),k) = NaN;
                numResponses(ceil(j/numPulses),k) = 0;
%                 spikeTimes_light_TT(ceil(j/numPulses),k) = NaN;
%                 spikeTimes_light_TT = [spikeTimes_light_TT tempSpike];
            end
            
            if k < numPulses/2 + 1
                laserSham = lL + (k-1)*(tB/(numPulses/2));
            end
            if k > numPulses/2
                laserSham = uL - (k-1)*(tB/(numPulses/2));
            end
            tempSpikeSham = laserStruct.(spikeFields{i})(find(laserStruct.(spikeFields{i}) > laserSham, 1));
            tempLatSham = tempSpikeSham - laserSham;
            if tempLatSham < 30 %30ms limit to evoke a spike
                latencyToResponseSham(ceil(j/numPulses),k) = tempLatSham - p.Results.TTdelay/1000;
                numResponsesSham(ceil(j/numPulses),k) = 1;
            else
                latencyToResponseSham(ceil(j/numPulses),k) = NaN;
                numResponsesSham(ceil(j/numPulses),k) = 0;
            end
        end
        if j == 1
            spikeTimes_spont_TT = [spikeTimes_spont_TT laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i})  < lL)];
        elseif j == trainInd(end)
            spikeTimes_spont_TT = [spikeTimes_spont_TT laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i})  >= uL)];
        else
            nuL = laserStruct.laserOff(j+numPulses) - tB;
            spikeTimes_spont_TT = [spikeTimes_spont_TT laserStruct.(spikeFields{i})(laserStruct.(spikeFields{i}) >= uL & laserStruct.(spikeFields{i}) < nuL)];
        end
    end
    spikeTimes_light_TT = 1000*spikeTimes_light_TT - p.Results.spikeDelay;
    
    if ~isempty(spikeTimes_light_TT)
        if p.Results.filtFlag == 1 
            if ~exist(strcat(filtDataPath,CSCchann{1},'.ncs'),'file')       %if filtered traces don't exit run function to filter CSC traces
                filterCSCtraces(unsortedDataPath, CSCchann, spikeTimes_light_TT);
                unsortedDataPath = filtDataPath;                            %change path to draw from filtered traces
            else
               unsortedDataPath = filtDataPath;
            end
        end
        
        [spikeTimes_light_CSC, spikeCSC_light] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{1},'.ncs'),[1 0 0 0 1],0,5,spikeTimes_light_TT);
        if length(spikeTimes_light_TT) > length(spikeTimes_light_CSC) % sometimes, spikeTimes_light_TT has indices that don't correspond to CSC data
            m = 1;
            while length(spikeTimes_light_TT) > length(spikeTimes_light_CSC)
                if spikeTimes_light_TT(m) - spikeTimes_light_CSC(m) < 0
                    spikeTimes_light_TT(m) = [];
                    m = m - 1;
                end
                m = m + 1;
            end
        end  
        
        [spikeCSC_light(:,:,2)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{2},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_light_TT);
        [spikeCSC_light(:,:,3)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{3},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_light_TT);
        [spikeCSC_light(:,:,4)] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{4},'.ncs'),[0 0 0 0 1],0,5,spikeTimes_light_TT);
        spikeInd_light = floor((spikeTimes_light_TT - spikeTimes_light_CSC)/31.25);
        
        [~, chanMaxInd] = max(max(abs([spikeCSC_light(:,1,1) spikeCSC_light(:,1,2) spikeCSC_light(:,1,3) spikeCSC_light(:,1,4)])));
        time = [(laserStruct.laserOn(1) - 200), (laserStruct.laserOn(numPulses) + 200)] * 1000;
        [timeStampTrace, exTrace] = Nlx2MatCSC(strcat(unsortedDataPath, CSCchann{chanMaxInd},'.ncs'),[1 0 0 0 1],0,4, time);
        timeDiff = (time(1) - timeStampTrace(1)) / 1000;
        lightTimes = [0:100:1000];
        lightTimes = lightTimes(1:end-1);
        lightEndTimes = lightTimes + 10;
        timePoints = [0:1000/32000:(length(exTrace(:))/32000)*1000] - 200 - timeDiff;
        timePoints = timePoints(1:end-1);
        figure; plot(timePoints, exTrace(:)*sCon, 'k'); hold on;
        for tmp = 1:length(lightTimes)
            plotShaded([lightTimes(tmp) lightEndTimes(tmp)], [-400 -400; 200 200], 'c')
        end
        set(gca, 'tickdir', 'out')
        ylabel('\muV')
        xlabel('Time (ms)')
        xlim([-200 1000])
        title([sessionFolder ' ' spikeFields{i}], 'Interpreter', 'none')
    end 
end
end