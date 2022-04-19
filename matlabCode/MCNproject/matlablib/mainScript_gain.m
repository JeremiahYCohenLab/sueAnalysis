% about FI curve, a, b slope and turning point; d how sharp the turning is 
I = 0.001*(1:2000);
F = FIcurve(I);
figure;
plot(I, F)
%% d
d = 0.05:0.02:0.2;
F = ones(size(d));
figure;hold on;
colors = cool(length(d));
for i = 1:length(d)
    F = FIcurve(I,'d',d(i));
    plot(I,F,'color',colors(i,:),'linewidth',1);
end
%% a
a = 200:5:300;
F = ones(size(a));
figure;hold on;
colors = cool(length(a));
for i = 1:length(a)
    F = FIcurve(I,'a',a(i));
    plot(I,F,'color',colors(i,:),'linewidth',1);
end

%% b
b = 200:5:300;
F = ones(size(b));
figure;hold on;
colors = cool(length(b));
for i = 1:length(b)
    F = FIcurve(I,'b',b(i));
    plot(I,F,'color',colors(i,:),'linewidth',1);
end
%% model parameters
modelparams = struct( ...
    'gE', 0.2609, ...
    'gI', -0.0497, ... % cross-inhibition strength [nA] 
    'I0', 0.3255, ... % background current [nA]
    'tauS', 0.1, ... % Synaptic time constant [sec]
    'gamma', 0.641, ... % Saturation factor for gating variable
    'tau0', 0.002, ... % Noise time constant [sec]
    'sigma', 0.02, ... % Noise magnitude [nA]
    'mu0', 20., ... % Stimulus firing rate [Hz]
    'Jext', 0.52, ... % Stimulus input strength [pA/Hz]
    'Ttotal', 2., ... % Total duration of simulation [s]
    'Tstim_on', 0.1, ... % Time of stimulus onset
    'Tstim_off', 1.6, ... % Time of stimulus offset
    'dt', 0.5/1000, ... % Simulation time step
    'record_dt', 5/1000, ...
    'decisionThresh', 8, ...
    'threshWin', 50/1000, ...% in s
    'a', 270); 
%% simulation with coherence
coh = 20;
s = simulation(modelparams,coh,5);
%% plotSingle trial
coh = 0;
plotLen = find(s.r1(:,1)~=0, 1, 'last' );
figure;
subplot(2,2,1); hold on;
plot(s.t(1:plotLen,:),s.I1(1:plotLen,:),'color','r','linewidth', 1.5)
plot(s.t(1:plotLen,:),s.I2(1:plotLen,:),'color','k', 'linewidth', 1.5)
title(['stimulus' ' coh' num2str(coh)])
subplot(2,2,2); hold on;
plot(s.t(1:plotLen,:),s.r1(1:plotLen,:),'color','r')
plot(s.t(1:plotLen,:),s.r2(1:plotLen,:),'color','k')
%line([modelparams.Tstim_on modelparams.Tstim_on], [0 1.1*max([s.r1; s.r2],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
%line([modelparams.Tstim_off modelparams.Tstim_off], [0 1.1*max([s.r1; s.r2],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
title('firing rate')
%set(gca, 'XScale', 'log')
subplot(2,2,3); hold on;
plot(s.t(1:plotLen),s.S1(1:plotLen,:),'color','r')
plot(s.t(1:plotLen),s.S2(1:plotLen,:),'color','k')
line([modelparams.Tstim_on modelparams.Tstim_on], [0 1.1*max([s.S1; s.S2],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
line([modelparams.Tstim_off modelparams.Tstim_off], [0 1.1*max([s.S1; s.S2],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
title('NMDA dynamics')
subplot(2,2,4); hold on;
plot(s.t(2:plotLen,:),s.Ieta1(2:plotLen,:),'color','r')
plot(s.t(2:plotLen,:),s.Ieta2(2:plotLen,:),'color','k')
line([modelparams.Tstim_on modelparams.Tstim_on], [min([s.Ieta1(2:end,:); s.Ieta2(2:end,:)],[],'all') max([s.Ieta1(2:end,:); s.Ieta2(2:end,:)],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
line([modelparams.Tstim_off modelparams.Tstim_off], [min([s.Ieta1(2:end,:); s.Ieta2(2:end,:)],[],'all') max([s.Ieta1(2:end,:); s.Ieta2(2:end,:)],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
title('noise')
%% Initialization
modelparams = struct( ...
    'gE', 0.08, ...
    'gI', -0.30, ... % cross-inhibition strength [nA] 
    'I0', 0.3255, ... % background current [nA]
    'tauS', 0.1, ... % Synaptic time constant [sec]
    'gamma', 0.641, ... % Saturation factor for gating variable
    'tau0', 0.002, ... % Noise time constant [sec]
    'sigma', 0.02, ... % Noise magnitude [nA]
    'mu0', 20., ... % Stimulus firing rate [Hz]
    'Jext', 0.52, ... % Stimulus input strength [pA/Hz]
    'Ttotal', 2., ... % Total duration of simulation [s]
    'Tstim_on', 0.1, ... % Time of stimulus onset
    'Tstim_off', 2, ... % Time of stimulus offset
    'dt', 0.5/1000, ... % Simulation time step
    'record_dt', 5/1000, ...
    'decisionThresh', 8, ...
    'threshWin', 50/1000, ...% in s
    'a', 270); 
%% 
Icommon = 0.15;
Idiff = 0.01;
s = simulationDMTwoSti(modelparams,0.03,Icommon+Idiff,Icommon-Idiff,8,false);
%% all tests
tau0Test = [0.002, 0.005, 0.01, 0.02];
aTest = [260:10:290];
Icommon = [0.1,0.2]; %Icommon = 0.05:0.05:0.2;  
Idiff = -0.0625:0.005:0.0625;
%Idiff = -0.0325:0.005:0.0325;
sim = 50;
savepath = '/Volumes/ZS256/MCNproject/simulations/';
%% tau and gain
for T = 4

    choices = zeros(length(aTest),length(Icommon),length(Idiff));
    RT = cell(length(aTest),length(Icommon),length(Idiff));
    tests = cell(length(aTest), length(Icommon), length(Idiff), sim);
    errors = -sim*ones(length(aTest),length(Icommon),length(Idiff));
    modelparams.tau0 = tau0Test(T);
    %% simulation choice and RTs
    for j = 1:length(aTest)
        modelparams.a = aTest(j);
        for i = 1:length(Icommon)
            for h = 1:length(Idiff)
                choices(j,i,h) = 0;
                RT{j,i,h} = [];
                for simNum = 1:sim
                    s = struct;
                    while sum(double(contains(fieldnames(s),'choice')))==0
                        s = simulationDMTwoSti(modelparams,0.03,Icommon(i)-Idiff(h),Icommon(i)+Idiff(h),8,true);
                        errors(j,i,h) = errors(j,i,h) + 1;
                    end
                    tests{j,i,h,simNum} = [s.r1,s.r2];
                    choices(j,i,h) = choices(j,i,h) + s.choice;
                    RT{j,i,h} = [RT{j,i,h} s.RT];
                end
            end
        end
    end
    
    save([savepath 'tau0=' num2str(tau0Test(T)) ' aTest' '.mat'], 'errors', 'choices', 'gItest', 'Icommon', 'Idiff', 'modelparams', 'RT', 'sim', 'aTest', 'tests', 'tau0Test')
end
%% simulate without choices
tests = cell(length(gItest), length(Icommon), length(Idiff));
for j = 1
    modelparams.gI = gItest(j);
    for i = 1:length(Icommon)
        for h = 1:length(Idiff)
                s = simulationDMTwoSti(modelparams,0.03,Icommon(i)-Idiff(h),Icommon(i)+Idiff(h),8,false);
                tests{j,i,h} = [s.r1,s.r2];
        end
    end
end
%% plotting all simulated traces
for j = 1
    figure; 
    for i = 1:length(Icommon)
        for h = 1:length(Idiff)
            frTemp = tests{j,i,h};
            subplot(length(Icommon),length(Idiff),(length(Idiff))*(i-1)+h); hold on;
            plot(s.t,frTemp(:,1),'color','r')
            plot(s.t,frTemp(:,2),'color','k')
            xlabel(['Idiff = ' num2str(Idiff(h))]);
            ylabel(['Icommon = ' num2str(Icommon(i))]);
        end
    end
    title(['gI=' num2str(gItest(j))])
end

%% plotting all choices and RTs
for h = 1:length(aTest)
figure;
imagesc(squeeze(choices(h,:,:)));
title(['a = ' num2str(aTest(h))])
figure;
subplot(1,2,1);
imagesc(cellfun(@mean,squeeze(RT(h,:,:))));
title(['a = ' num2str(aTest(h))])
subplot(1,2,2);
imagesc(cellfun(@std,squeeze(RT(h,:,:))));
title(['a = ' num2str(aTest(h))])
end
%% plot action selection sep by Icommon & Idiff

for h = 1:length(aTest)
    figure; hold on;
    colors = cool(size(choices,2));
    currChoices = squeeze(choices(h,:,:));
    for i = 1:size(choices,2)
        p = plot(Idiff, 0.5*(currChoices(i,:)/sim + 1), 'linewidth', 1.5, 'color', colors(i,:));
    end
    legend(cellfun(@num2str,num2cell(Icommon), 'UniformOutput', false))
    title(['p(choice == 1)' 'a = ' num2str(aTest(h))]);
end
%% plot action selection sep by Idiff & gain
colors = cool(size(choices,1));
figure; hold on;
for h = 1:length(aTest)
    currChoices = squeeze(choices(h,:,:));
    currChoices = mean(currChoices,1);
    plot(Idiff, 0.5*(currChoices/sim + 1), 'linewidth', 1.5, 'color', colors(h,:));
end
legend(cellfun(@num2str,num2cell(aTest), 'UniformOutput', false))
%%
figure;
title('p(choice == 1)'); hold on;
colors = cool(size(choices,1));
for i = 1:size(choices,1)
    plot(2*Idiff, 0.5*(choices(i,:)/20 + 1), 'linewidth', 1.5, 'color', colors(i,:));
end
legend({num2str(Isum(1))})
%%
figure; hold on;
for i = 1:size(RT,1)
    for j = 1:size(RT,2)
        subplot(size(RT,1), size(RT,2), size(RT,2)*(i-1)+j);
        histogram(RT{i,j},'facecolor',colors(i,:));       
    end
end
%%
gEtest = 0.01:0.05:0.9;
gItest = -0.5:0.05:-0.02;
%%
I0test = 0.3255;
figure; hold on;
    figure; hold on;
    for i = 1:length(gEtest)
        modelparams.gEtest = gEtest(i);
        for j = 1:length(gItest)
            subplot(length(gEtest), length(gItest), length(gItest)*(i-1)+j);hold on;
            modelparams.gI = gItest(j);
            s = simulationDM(modelparams,0.33,0.27,8,false);
            plot(s.t,s.r1,'color','r')
            plot(s.t,s.r2,'color','k')
        end
    end
%%
% inputs 0.33 0.27
% inhibition 0.05
% I0 cannot be too big
% worked: s = simulationDM(modelparams,0.1,0.08,0.05,8,false); gE 0.15 gI -0.3


%% fit all choice curves for gain
tau0Test = [0.02  0.01, 0.005 0.002];
sigma = zeros(length(tau0Test), length(aTest),3);
for i = 1:length(tau0Test)
    load([savepath 'tau0 =' num2str(tau0Test(i)) ' aTest' '.mat'])
    x = Idiff*2;
     for j = 1:length(aTest)
        currChoices = squeeze(choices(j,:,:));
        currChoices = mean(currChoices,1);
        y = 0.5*(currChoices/sim + 1);
        currFit = fitSigmoid(x',y');
        sigma(i,j,1) = currFit.a;
        sigma(i,j,2:3) = confint(currFit);
     end
end
figure;
imagesc(sigma(:,:,1));
%%
figure; hold on;
colors = cool(length(tau0Test));
for i = 1:length(tau0Test)
    errorbar(gItest,squeeze(sigma(i,:,1)),0.5*(squeeze(sigma(i,:,3)-sigma(i,:,2))), 'LineWidth', 2, 'Color', colors(i,:))
end
legend(cellfun(@num2str,num2cell(tau0Test), 'UniformOutput', false));
%%
figure; hold on;
colors = cool(length(gItest));
for i = 1:length(gItest)
    errorbar(tau0Test,squeeze(sigma(:,i,1)),0.5*(squeeze(sigma(:,i,3)-sigma(:,i,2))), 'LineWidth', 2, 'Color', colors(i,:))
end
legend(cellfun(@num2str,num2cell(gItest), 'UniformOutput', false));
%% choose tau0=0.01 as standard
i = find(tau0Test == 0.01);
lm = fitlm(gItest',sigma(i,:,1)');
lm.Coefficients.Estimate
%%