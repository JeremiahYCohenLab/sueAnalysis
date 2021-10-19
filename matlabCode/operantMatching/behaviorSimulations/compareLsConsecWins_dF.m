function compareLsConsecWins_dF(xlFile, sheet, category, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('modelNames', {'fiveParam_bias', 'eightParam_absPePeAN_scale_int_bias'});
p.addParameter('beh', 'clean')
p.addParameter('bernFlags', [])
p.addParameter('runs', 500)
p.addParameter('rSeed', 56715)
p.addParameter('rwdProbs', [90 50 10])
p.parse(varargin{:});

if isempty(p.Results.bernFlags)
    bernFlags = ones(1,length(p.Results.modelNames));
else
    bernFlags = p.Results.bernFlags;
end

rSeed = p.Results.rSeed;

[root, sep] = currComputer();

numMdls = length(p.Results.modelNames);

[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end


prevAnimal = [];
rwds = [];
changeChoice = [];
choices = [];
lS = [];
for currS = 1:length(dayList)
    sessionName = dayList{currS};
    [animal, date] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    date = date(1:9);
    sessionFolder = ['m' animal date];

    if strcmp(animal, prevAnimal) == 0
      %  runFrac = sum(~cellfun(@isempty,strfind(dayList, animal)) == 1) / length(dayList);
      %  runs = ceil(p.Results.runs * runFrac);
        if isempty(prevAnimal)
            lS_sim = nan(numMdls,10);
            a = 1;
        else
            lS_sim = cat(3, lS_sim, nan(numMdls,10));
            a = a+1;
        end
        for currM = 1:numMdls
            if bernFlags(currM)
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNames{currM}...
                    sep animal p.Results.beh '_' p.Results.modelNames{currM} '.mat'];
            else
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelNames{currM}... 
                    sep animal p.Results.beh '_' p.Results.modelNames{currM} '.mat'];
            end
            t = generateStanModelTerms_opMD(p.Results.modelNames{currM}, modelPath, [], 0);
            params = t.params;

            simRwds = []; simChangeChoice = []; simChoices = [];
            for currR = 1:p.Results.runs
                [~, rwdTmp, choiceTmp, ~, ~] = runSim_dF(p.Results.modelNames{currM}, params,...
                            300, rSeed, p.Results.rwdProbs);
                rSeed = rSeed + 1;
                simRwds = [simRwds nan(1,10) abs(rwdTmp)];
                simChangeChoice = [simChangeChoice nan(1,10) (abs(diff(choiceTmp)) > 0) NaN];
                simChoices = [simChoices nan(1,10) choiceTmp];
            end

            for i = 1:10
                tmpPat = [0 ones(1,i) 0];               %find instances of outcome pattern
                inds = strfind(simRwds, tmpPat);        
                tmpPatC = ones(1,length(tmpPat)-1);
                indsC = strfind(simChoices, tmpPatC);   %find instances of consecutive choices
                tmpPatCn = ones(1,length(tmpPat)-1) * -1;
                indsCn = strfind(simChoices, tmpPatCn);
                indsC = sort(unique([indsC indsCn])) - 1;
                
                if ~isempty(inds)
                    inds = inds(ismember(inds, indsC));                         %limit outcome pattern to instances of consecutive choices
                    inds(isnan(simChangeChoice(inds+length(tmpPat)-1))) = [];   %get rid of choices at end of session
                    lS_sim(currM, i, a) = nansum(simChangeChoice(inds+length(tmpPat)-1)) / length(inds);
                end
            end
        end
        prevAnimal = animal;
        newAnimal = 1;
    end

    if (newAnimal & a > 1) || currS == length(dayList)
        lS = [lS; nan(1,10)];
        for i = 1:10
            tmpPat = [0 ones(1,i) 0];               %find instances of outcome pattern
            inds = strfind(rwds, tmpPat);
            tmpPatC = ones(1,length(tmpPat)-1);
            indsC = strfind(choices, tmpPatC);      %find instances of consecutive choices
            tmpPatCn = ones(1,length(tmpPat)-1) * -1;
            indsCn = strfind(choices, tmpPatCn);
            indsC = sort(unique([indsC indsCn])) - 1;

            if ~isempty(inds)
                inds = inds(ismember(inds, indsC));                         %limit outcome pattern to instances of consecutive choices
                inds(isnan(changeChoice(inds+length(tmpPat)-1))) = [];      %get rid of choices at end of session
                lS(a-1, i) = nansum(changeChoice(inds+length(tmpPat)-1)) / length(inds);
            end
        end
        rwds = []; changeChoice = []; choices = [];
        newAnimal = 0;
    end
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData([sessionName '.asc'], 0);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
    rwds = [rwds nan(1,10) abs(behavStruct.allRewards)];
    changeChoice = [changeChoice nan(1,10) (abs(diff(behavStruct.allChoices)) > 0) NaN];
    choices = [choices nan(1,10) behavStruct.allChoices];
end
    


figure; hold on;
colors = cool(numMdls + 1);
x = [1:10];
plotFilled(x, lS, colors(1,:));
legTxt = [];
for currM = 1:numMdls
    plotFilled(x, squeeze(lS_sim(currM,:,:))', colors(currM+1,:))
    legTxt = [legTxt {strrep(p.Results.modelNames{currM}, '_', ' ')} {''}];
end
xlabel('consecutive wins before loss')
ylabel('lose-shift')
legend([{'actual'} {''} legTxt]) 

set(gca, 'tickdir' ,'out')
set(gcf, 'renderer', 'painters')

end
