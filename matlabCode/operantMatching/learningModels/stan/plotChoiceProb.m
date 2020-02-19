function [model] = plotChoiceProb(sessionName, beh, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('modelNames', {'sixParam_absPePeAN_bi'})
p.addParameter('plotFlag', 1);
p.addParameter('bernFlag', [1]);
p.addParameter('sessionFlag', 0);
p.parse(varargin{:});

[root, sep] = currComputer();
modelNames = p.Results.modelNames;

filename = [sessionName '.asc'];
[behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';
ITI = [behavStruct.timeBtwn]';


if p.Results.plotFlag
    normKern = normpdf(-15:15,0,4);
    normKern = normKern / sum(normKern);
    normKern = fliplr(normKern(1:16));
    allChoices = choice(:,1) - choice(:,2);
    figure; hold on;
    yyaxis left; plot(conv(allChoices',normKern)/max(conv(allChoices',normKern)),'k','linewidth',1);
    ylim([-1 1]); yticks([-1 0 1]); yticklabels([0 0.5 1])
    set(gca, 'Ycolor', 'k');
    
    numParams = length(modelNames);
    blue = [0 1 1];
    purp = [0.7 0 1];
    colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
end


[animalName, date] = strtok(sessionName, 'd'); 
animalName = animalName(2:end);

for currMod = 1:length(modelNames)
    
    if p.Results.bernFlag(currMod)
        modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelNames{currMod} sep animalName...
        beh '_' modelNames{currMod} '.mat'];
    else
        modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep modelNames{currMod} sep animalName...
        beh '_' modelNames{currMod} '.mat'];
    end
    
    t = generateStanModelTerms_opMD(modelNames{currMod}, modelPath, sessionName, p.Results.sessionFlag, p.Results.revForFlag);
    
    model.(modelNames{currMod}) = t;
    
    if p.Results.plotFlag
        yyaxis right; 
        eval(['plot(model.' modelNames{currMod} '.probChoice(:,1), ''-'', ''Color'', [' num2str(colors(currMod,:)) '],''linewidth'', 1)']);
        ylim([0 1]); yticks([0 0.5 1]);
        set(gca, 'Ycolor', 'k');
    end 
    
end

if p.Results.plotFlag
    ax = gca;
    ax.TickDir = 'out';
    set(gcf, 'Position', [-1910 159 1905 750]);
    set(gcf, 'Renderer', 'Painters')
    legend(['choices', modelNames], 'Interpreter', 'none')
    suptitle(sessionName)
end



