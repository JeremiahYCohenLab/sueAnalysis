function [t, o, noSession] = getStanModelParams_samps(modelName, modelPath, numSamps, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('biasFlag',1)
p.addParameter('sessionParamsFlag', 1)
p.addParameter('sessionName', [])
p.addParameter('varFlag', 1)
p.addParameter('revForFlag', 0)
p.parse(varargin{:});

paramNames = getParamNames_dF(modelName, p.Results.biasFlag);

%get model params
load(modelPath);
if p.Results.sessionParamsFlag
    sessionInd = find(~cellfun(@isempty,strfind(dayList,p.Results.sessionName)));
end
noSession = false;
t = struct;
o = struct;
if isempty(sessionInd)
    noSession = true;
    return
end
tmp = whos;
samples = eval(tmp(1).name);

params = nan(numSamps, length(paramNames));
inds = randperm(length(samples.lp__));
inds = inds(1:numSamps);
for currP = 1:length(paramNames)
    if p.Results.sessionParamsFlag
        tmp = eval(['samples.' paramNames{currP}]);
        tmp = tmp(:,sessionInd);
    else
        if ~strcmp(paramNames{currP}, 'bias')
            tmp = eval(['samples.mu_' paramNames{currP}]);
        end
    end
    params(:,currP) = tmp(inds);
end


t = struct;
t.params = params;
if p.Results.varFlag
    %get session behavior
    [behSessionData,blockSwitch,~] = loadBehavioralData([p.Results.sessionName '.asc'], p.Results.revForFlag);
    o = parseBehavioralData(behSessionData, blockSwitch);
    outcome = abs([o.allRewards])';
    choice = o.allChoices';
    choice(o.allChoices<0) = 0;
    
    for currS = 1:numSamps
        tmpStruct{currS} = getModelVariables_dF(modelName, params(currS, :), choice, outcome);
        
        %get rid of any bad parameters that cause outlier variable values
        while any(abs(tmpStruct{currS}.pe) > 1)
            tmpInd = randperm(length(samples.lp__));
            while any(inds == tmpInd(1))
                tmpInd = randperm(length(samples.lp__));
            end
            for currP = 1:length(paramNames)
                if p.Results.sessionParamsFlag
                    tmp = eval(['samples.' paramNames{currP}]);
                    tmp = tmp(:,sessionInd);
                else
                    tmp = eval(['samples.mu_' paramNames{currP}]);
                end
                params(currS,currP) = tmp(tmpInd(1));
            end 
            tmpStruct{currS} = getModelVariables_dF(modelName, params(currS, :), choice, outcome);
        end
    end
    
    infInds = [];
    mdlVarNames = fields(tmpStruct{1});
    for currV = 1:length(mdlVarNames)
        if isempty(regexp(mdlVarNames{currV}, 'Q'))
            tmp = [];
            for currS = 1:numSamps
                tmp = [tmp tmpStruct{currS}.(mdlVarNames{currV})];
            end
            if sum(sum(isinf(tmp))) > 0
                [~,c] = find(isinf(tmp));
                infInds = [infInds, c'];
            end
            tmp(:,infInds) = NaN;
            t.(mdlVarNames{currV}) = nanmean(tmp,2); 
        end
    end
 
    tmpQ_L = []; tmpQ_R =[];
    tmpP_L = []; tmpP_R =[];
    for currS = 1:numSamps
        tmpQ_L = [tmpQ_L tmpStruct{currS}.Q(:,1)];
        tmpQ_R = [tmpQ_R tmpStruct{currS}.Q(:,2)];
    end
    tmpQ_L(:,infInds) = NaN;
    tmpQ_R(:,infInds) = NaN;
    t.Q(:,1) = nanmean(tmpQ_L,2);
    t.Q(:,2) = nanmean(tmpQ_R,2);
else
    o = [];
end
    
    
    
    
    

