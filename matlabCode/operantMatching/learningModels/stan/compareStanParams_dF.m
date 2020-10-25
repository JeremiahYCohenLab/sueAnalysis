function compareStanParams_dF(animals, pre, post, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('paramNames', {'aNscale', 'aNmin', 'aP', 'aF', 'aPE', 'beta'});
p.addParameter('modelName', ['sixParam_absPePeAN_bi']);
p.addParameter('hyperFlag', 1);
p.addParameter('bernFlag', 1);
p.parse(varargin{:});

[root, sep] = currComputer();
numParams = length(p.Results.paramNames);
numAnimals = length(animals);

for aI = 1:length(animals)
    if p.Results.hyperFlag
        if p.Results.bernFlag
            filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep p.Results.modelName '_manip' sep];
        else
            filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName '_manip' sep];
        end
        mdlFolder = dir(filePath);
        fileInd = (~cellfun(@isempty,strfind({mdlFolder.name}, [p.Results.modelName '_manip.mat'])));
        if sum(fileInd) ~= 1
            uiopen(filePath)
            tmp = whos('-regexp', animals{aI});
            samps = eval(tmp.name);
        else
            load([filePath mdlFolder(fileInd).name]);
            tmp = whos('-regexp', animals{aI});
            samps = eval(tmp.name);
        end
        
        for ind = 1:length(p.Results.paramNames);
            avgParamsPre{ind}(aI, 1) = eval(['mean(samps.mu_' p.Results.paramNames{ind} ')']);
            avgParamsPre{ind}(aI, 2) = eval(['std(samps.mu_' p.Results.paramNames{ind} ')']);
            avgParamsPost{ind}(aI, 1) = eval(['mean(samps.d_mu_' p.Results.paramNames{ind} ')']);
            avgParamsPost{ind}(aI, 2) = eval(['std(samps.d_mu_' p.Results.paramNames{ind} ')']); 
        end

    else
        filePath = [root sep animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep p.Results.modelName sep];
        mdlFolder = dir(filePath);
        tmpInd = (~cellfun(@isempty,strfind({mdlFolder.name}, pre)));
        load([filePath mdlFolder(tmpInd).name]);
        sampsPre = eval([animals{aI} pre '_' p.Results.modelName]);
        tmpInd = (~cellfun(@isempty,strfind({mdlFolder.name}, post)));
        load([filePath mdlFolder(tmpInd).name]);
        sampsPost = eval([animals{aI} post '_' p.Results.modelName]);
        
        for ind = 1:numParams
            avgParamsPre{ind}(aI, 1) = eval(['mean(sampsPre.mu_' p.Results.paramNames{ind} ')']);
            avgParamsPre{ind}(aI, 2) = eval(['std(sampsPre.mu_' p.Results.paramNames{ind} ')']);
            avgParamsPost{ind}(aI, 1) = eval(['mean(sampsPost.mu_' p.Results.paramNames{ind} ')']);
            avgParamsPost{ind}(aI, 2) = eval(['std(sampsPost.mu_' p.Results.paramNames{ind} ')']); 
        end
    end
    
end

blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numAnimals)', linspace(blue(2),purp(2),numAnimals)', linspace(blue(3),purp(3),numAnimals)'];

titles = generateParamTitles(p.Results.paramNames);

figure;
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    for j = 1:length(animals)
        errorbar([avgParamsPre{i}(j,1) avgParamsPost{i}(j,1)], [avgParamsPre{i}(j,2) avgParamsPost{i}(j,2)],...
            'Color', colors(j,:), 'linewidth', 1.5);
    end
    plot([mean(avgParamsPre{i}(:,1)) mean(avgParamsPost{i}(:,1))], 'Color', 'k', 'linewidth', 4)
    xticks([1 2])
    xticklabels({'pre', 'post'})
    xlim([0.5 2.5])
    set(gca,'tickdir', 'out')
    if isempty([strfind(p.Results.paramNames{i}, 'beta') strfind(p.Results.paramNames{i}, 'k') strfind(p.Results.paramNames{i}, 'bias')])
        ylim([0 1])
    end
    title(titles{i})
end
legend(animals)
titleTxt = strrep([p.Results.modelName], '_', ' ');
if p.Results.hyperFlag
    titleTxt = [titleTxt ' hyperhyper'];
else
    titleTxt = [titleTxt ' separate fits'];
end
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')
