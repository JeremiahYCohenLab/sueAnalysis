function compareDynaVolTransition(xlFile, sheet, col, varargin)
    p = inputParser;
    % default parameters if none given
    p.addParameter('revForFlag', 0);
    p.addParameter('numTrialPre', 5);
    p.addParameter('numTrialPost', 10);
    p.addParameter('halves', 0);
    p.parse(varargin{:});
    
    switch p.Results.halves 
        case 0
            position = 'all';
        case 1
            position = 'early';
        case 2
            position = 'late';
    end


    
    dayList = getDayList(xlFile, sheet, col);
    
    choicesAroundTrans = [];
    choicesAtJunction = [];
    combinedVol = [];
    combinedHalves = [];
    combinedFirstHalf = [];
    
 

for i = 1:length(dayList)
    sessionName = dayList{i};
    s = behAnalysisNoPlot_opMD(sessionName, 'simpleFlag', 1);
    volChange = find(s.vol(1:end-1)~=s.vol(2:end))+1;
    allChoices = s.allChoices;
    allChoices(allChoices<0) = 0;
    for j = 2:length(s.blockSwitch)
        startTrial = max(0, s.blockSwitch(j)-p.Results.numTrialPre);
        endTrial = min(length(s.allChoices), s.blockSwitch(j)+p.Results.numTrialPost-1);
        tmpChoices = [NaN(1,-s.blockSwitch(j)+startTrial+p.Results.numTrialPre) ...
                      allChoices(startTrial:endTrial) ...
                      NaN(1,s.blockSwitch(j)+p.Results.numTrialPost-1-endTrial)];
        if s.blockProbs(s.blockSwitch(j),1) < 50
            tmpChoices = 1-tmpChoices;
        end
        if s.vol(s.blockSwitch(j)) == s.vol(s.blockSwitch(j)-1) %% discard those at transitions
            choicesAroundTrans = [choicesAroundTrans; tmpChoices];
            combinedVol = [combinedVol; s.vol(s.blockSwitch(j))]; 
            combinedHalves = [combinedHalves; double(s.blockSwitch(j)>volChange)];
            combinedFirstHalf = [combinedFirstHalf; s.vol(1)];
        else
            choicesAtJunction = [choicesAtJunction; tmpChoices];
        end
    end
end
    figure;
    subplot(3,4,1); hold on;
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    title([sheet ' ' col 'all'])
    subplot(3,4,2); hold on;
    allSameChoiceInds = sum(choicesAroundTrans(:,1:p.Results.numTrialPre+1), 2);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&allSameChoiceInds==p.Results.numTrialPre+1,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&allSameChoiceInds==p.Results.numTrialPre+1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    
    
    subplot(3,4,3); hold on;
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&combinedHalves==1,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&combinedHalves==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    title([sheet ' ' col 'Second'])
    subplot(3,4,4); hold on;
    allSameChoiceInds = sum(choicesAroundTrans(:,1:p.Results.numTrialPre+1), 2);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&allSameChoiceInds==p.Results.numTrialPre+1&combinedHalves==1,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&allSameChoiceInds==p.Results.numTrialPre+1&combinedHalves==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    

    subplot(3,4,5); hold on;
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&combinedHalves==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&combinedHalves==0,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    title([sheet ' ' col 'First'])
    subplot(3,4,6); hold on;
    allSameChoiceInds = sum(choicesAroundTrans(:,1:p.Results.numTrialPre+1), 2);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&allSameChoiceInds==p.Results.numTrialPre+1&combinedHalves==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&allSameChoiceInds==p.Results.numTrialPre+1&combinedHalves==0,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    
    
    subplot(3,4,7); hold on;
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&combinedFirstHalf==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&combinedFirstHalf==0,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'})
    title([sheet ' ' col  'startWithLow']);
    subplot(3,4,8); hold on;
    allSameChoiceInds = sum(choicesAroundTrans(:,1:p.Results.numTrialPre+1), 2);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&allSameChoiceInds==p.Results.numTrialPre+1&combinedFirstHalf==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&allSameChoiceInds==p.Results.numTrialPre+1&combinedFirstHalf==0,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    
  
    subplot(3,4,9); hold on;
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&combinedFirstHalf==1,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&combinedFirstHalf==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    title([sheet ' ' col 'startWithHigh'])
    subplot(3,4,10); hold on;
    allSameChoiceInds = sum(choicesAroundTrans(:,1:p.Results.numTrialPre+1), 2);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==0&allSameChoiceInds==p.Results.numTrialPre+1&combinedFirstHalf==1,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedVol==1&allSameChoiceInds==p.Results.numTrialPre+1&combinedFirstHalf==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    
    
    subplot(3,4,11); hold on;
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedFirstHalf==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(combinedFirstHalf==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    title([sheet ' ' col 'startWith'])
    subplot(3,4,12); hold on;
    allSameChoiceInds = sum(choicesAroundTrans(:,1:p.Results.numTrialPre+1), 2);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(allSameChoiceInds==p.Results.numTrialPre+1&combinedFirstHalf==0,:), [0 0 1]);
    plotFilledBern(-p.Results.numTrialPre:1:p.Results.numTrialPost-1, choicesAroundTrans(allSameChoiceInds==p.Results.numTrialPre+1&combinedFirstHalf==1,:), [1 0 1]);
    legend({'low volatility', '', 'high volatlity'});
    
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
    
    
    
    

    

    
    
