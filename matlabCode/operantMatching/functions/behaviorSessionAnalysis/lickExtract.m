function struct = lickExtract(session, plotFlag)
os = behAnalysisNoPlot_opMD(session, 'simpleFlag',1);
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
load([pd.sortedFolder session '_sessionData_behav.mat']);
solenoidTime = 1800 + os.rwdDelay + 1000; % response window, extended rwd time, preITI time
licksTrialL = cell(length(behSessionData),1);
licksTrialR = cell(length(behSessionData),1);
licksTrialAll = cell(length(behSessionData),1);
licksITIL = cell(length(behSessionData),1);
licksITIR = cell(length(behSessionData),1);
licksITIAll = cell(length(behSessionData),1);
side = zeros(length(behSessionData),1);
sideITI = zeros(length(behSessionData),1);
itiLickTimeL = NaN(length(behSessionData),1);
itiLickTimeR = NaN(length(behSessionData),1);
itiLickTimeAll = NaN(length(behSessionData),1);
for i = 1:length(behSessionData)
    tmpL = behSessionData(i).licksL - behSessionData(i).CSon;
    tmpL = tmpL(tmpL <= solenoidTime);
    tmpR = behSessionData(i).licksR - behSessionData(i).CSon;
    tmpR = tmpR(tmpR <= solenoidTime);
    licksTrialL{i} = tmpL;
    licksTrialR{i} = tmpR;
    licksTrialAll{i} = sort([tmpL tmpR]);

    if ~isempty(tmpL) || ~isempty(tmpR)
        if isempty(tmpL)
            side(i) = 1;
        else
            if isempty(tmpR)
                side(i)= -1;
            else
                [~, side(i)] = min([min(tmpL) min(tmpR)]);
                side(i) = 2*side(i)-3;
            end
        end
    end
    
    tmpL = behSessionData(i).licksL - behSessionData(i).CSon - solenoidTime;
    tmpL = tmpL(tmpL > 0);
    tmpR = behSessionData(i).licksR - behSessionData(i).CSon - solenoidTime;
    tmpR = tmpR(tmpR > 0); 
    licksITIL{i} = tmpL;
    licksITIR{i} = tmpR;
    licksITIAll{i} = sort([tmpL tmpR]);
    if ~isempty(tmpL) || ~isempty(tmpR)
        if isempty(tmpL)
            sideITI(i) = 1;
            itiLickTimeR(i) = min(tmpR);
            itiLickTimeAll(i) = min(tmpR);
        else
            if isempty(tmpR)
                sideITI(i) = -1;
                itiLickTimeL(i) = min(tmpL);
                itiLickTimeAll(i) = min(tmpL);
            else
                [~, sideITI(i)] = min([min(tmpL) min(tmpR)]);
                sideITI(i) = 3;
                itiLickTimeL(i) = min(tmpL);
                itiLickTimeR(i) = min(tmpR);
                itiLickTimeAll(i) = min([tmpL tmpR]);
            end
        end
    end
    
end
struct.solenoidTime = solenoidTime;
struct.side = side;
struct.sideITI = sideITI;

struct.licksITIL = licksITIL;
struct.licksITIR = licksITIR;
struct.licksITIAll = licksITIAll;

struct.licksTrialL = licksTrialL;
struct.licksTrialR = licksTrialR;
struct.licksTrialAll = licksTrialAll;

struct.itiLatL = itiLickTimeL;
struct.itiLatR = itiLickTimeR;
struct.itiLatAll = itiLickTimeAll;

sideITI(sideITI==3) = 0;
if plotFlag
    figure;
    subplot(3,1,1); hold on;
    plot([os.responseInds; os.responseInds], [zeros(length(os.allChoices),1), 0.5*os.allChoices' + 0.5*os.allRewards']', 'k');
    plot([os.responseInds+0.5; os.responseInds+0.5], [zeros(1,length(os.allChoices)); sideITI(os.responseInds)'], 'r', 'LineWidth',3);
    
    subplot(3,3,4);
    edges = linspace(min(itiLickTimeAll-1), max(itiLickTimeAll)+1, 20);
    histogram(itiLickTimeAll, edges, 'FaceColor', [0.6 0.6 0.6], 'Normalization', 'probability');
    title('ITI lick time in all trials')
    xlabel('time (ms)')
%     set(gca,'xscale','log')
    
    subplot(3,3,5);
    hold on;
    histogram(itiLickTimeL, edges, 'FaceColor', 'm', 'Normalization', 'probability');
    histogram(itiLickTimeR, edges, 'FaceColor', 'c', 'Normalization', 'probability');
    legend({'L', 'R'})
    title('ITI lick time on two sides')
%     set(gca,'xscale','log')
    xlabel('time (ms)')
    
    subplot(3,3,7);
    hold on;
    numITIAll = cellfun(@length, licksITIAll);
    numITIL = cellfun(@length, licksITIL);
    numITIR = cellfun(@length, licksITIR);
    edges = linspace(0.5,max(numITIAll)+0.01, 20);
    histogram(numITIAll, edges, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability')    
    title('ITI lick num in all trials')
%     set(gca,'xscale','log')
    xlabel('lick num')
    
    subplot(3,3,8);
    hold on;
    histogram(numITIL, edges, 'FaceColor', 'm', 'Normalization', 'probability')
    histogram(numITIR, edges, 'FaceColor', 'c', 'Normalization', 'probability')
    title('ITI lick num on two sides')
%     set(gca,'xscale','log')
    xlabel('lick num')
     
    subplot(3,3,9);
    csOn = [behSessionData.CSon];
    itiLen = [csOn(2:end) - csOn(1:end-1), NaN] -solenoidTime-1000;
    scatter(itiLen, numITIAll);
%     set(gca,'xscale','log')
    xlabel('ITI length ?ms?');
    ylabel('num licks')
    sgtitle(session);

    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
    
    saveFigurePDF(gcf,[pd.saveFigFolder session '_itiLick.pdf'])

end

end