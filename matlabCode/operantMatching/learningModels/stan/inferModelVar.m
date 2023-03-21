function t = inferModelVar(session, params, modelName, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('biasFlag',1)
p.addParameter('revForFlag', 0)
p.addParameter('perturb', []);
p.parse(varargin{:});

t = struct;
t.params = params;

%get session behavior
o = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
outcome = abs([o.allRewards])';
choice = o.allChoices';
choice(o.allChoices<0) = 0;
ITI = o.timeBtwn;
currInd = 1;
for currS = 1:size(params,1)
    if isempty(p.Results.perturb)
        tmp = getModelVariables_dF(modelName, params(currS, :), choice, outcome);
    else
        tmp = getModelVariablesLaser_dF(modelName, params(currS, :), choice, outcome, o.laser);
    end

        tmpStruct{currInd} = tmp; 
        currInd = currInd + 1;
    
end

infInds = [];
mdlVarNames = fields(tmpStruct{1});
for currV = 1:length(mdlVarNames)
    if isempty(regexp(mdlVarNames{currV}, 'Q', 'once'))
        tmp = [];
        for currS = 1:length(tmpStruct)
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
for currS = 1:length(tmpStruct)
    tmpQ_L = [tmpQ_L tmpStruct{currS}.Q(:,1)];
    tmpQ_R = [tmpQ_R tmpStruct{currS}.Q(:,2)];
end
tmpQ_L(:,infInds) = NaN;
tmpQ_R(:,infInds) = NaN;
t.Q(:,1) = nanmean(tmpQ_L,2);
t.Q(:,2) = nanmean(tmpQ_R,2);
    
    
    
    
    

