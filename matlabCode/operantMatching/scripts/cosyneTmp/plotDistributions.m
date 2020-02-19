%% plot manipulation parameters
paramNames = {'aNscale', 'aNmin', 'aP','aF', 'aPE', 'beta'};
paramTitles= generateParamTitles(paramNames);
numParams = length(paramNames);


figure; 
blue = [0 0.7 1];
purp = [0.5 0 0.8];
colors = [linspace(blue(1),purp(1),numParams*2)', linspace(blue(2),purp(2),numParams*2)', linspace(blue(3),purp(3),numParams*2)'];
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['meow.mu_' paramNames{i}]) , 50,...
        'Normalization', 'Probability', 'FaceColor', blue, 'EdgeColor', 'none')
    histogram(eval(['meow.mu_d_' paramNames{i}]) , 50,...
        'Normalization', 'Probability', 'FaceColor', purp, 'EdgeColor', 'none')
    set(gca,'tickdir', 'out') 
    title(paramTitles{i})
end
suptitle('CG84');
legend('pre', 'post')
set(gcf,'Renderer', 'Painters')


%% generate best estimates of parameters
paramEsts = nan(1,length(numParams));
for i = 1:numParams
    paramEsts(i)  = nanmedian(eval(['CG79clean_sixParam_absPePeAN.mu_' paramNames{i}]));
end

%% plot session parameters
numSesh = 6;
figure; 
suptitle('CG14 beta min');
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numSesh)', linspace(blue(2),purp(2),numSesh)', linspace(blue(3),purp(3),numSesh)'];
for i = 1:numSesh
    subplot(1,numSesh,i); hold on;
    histogram(samples_60.mu_(:,i), 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    set(gca,'tickdir', 'out') 
    title(['sesh ' num2str(i)])
end
set(gcf,'Renderer', 'Painters')


%% plot hyperparameters
paramNames = {'aNscale', 'aNmin', 'aP','aF', 'aPE', 'beta'};
numParams = length(paramNames);

figure; 
suptitle('CG75');
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['CG75clean_sixParam_absPePeAN_bi.mu_' paramNames{i}]) , 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    title(paramNames{i})
end
set(gcf,'Renderer', 'Painters')


%%
paramNames = {'aN', 'aP','aF', 'beta'};
numParams = length(paramNames);
numSessions = size(eval(['meow.' paramNames{1}]),2);
paramTitles = generateParamTitles(paramNames);

figure; 
suptitle('CG40 - fourParam');
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numSessions)', linspace(blue(2),purp(2),numSessions)', linspace(blue(3),purp(3),numSessions)'];
for i = 1:numParams
    subplot(3,numParams,[i i+numParams]); hold on;
    title(paramTitles{i})
    for j = 1:numSessions
        histogram(eval(['meow.' paramNames{i} '(:,j)']) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(j,:), 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    end
    set(gca, 'tickdir', 'out')
    subplot(3,numParams,(i+2*numParams)); hold on;
    paramTmp =[];
    for j = 1:numSessions
        tmp = eval(['meow.' paramNames{i} '(:,j)']);
        paramTmp = [paramTmp tmp];
    end
    plotFilledStd([1:j], paramTmp, [0 0 0]);
    set(gca, 'tickdir', 'out')
end
set(gcf,'Renderer', 'Painters')



