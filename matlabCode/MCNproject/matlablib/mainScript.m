% about FI curve, a, b slope and turning point; d how sharp the turning is 
I = 0.001*(1:1000);
F = FIcurve(I);
figure;
plot(I, F, 'linewidth',5, 'color', '.k')
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
    'Tstim_off', 1.0, ... % Time of stimulus offset
    'dt', 0.5/1000, ... % Simulation time step
    'record_dt', 5/1000, ...
    'decisionThresh', 8, ...
    'threshWin', 50/1000, ...% in s
    'a', 270); 
%% simulation with coherence
coh = 20;
s = simulation(modelparams,coh,5);
%% plotSingle trial
coh =0;
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
%%
figure; hold on;
plot(s.t(1:plotLen,:),s.r1(1:plotLen,:),'color','r','linewidth',3)
plot(s.t(1:plotLen,:),s.r2(1:plotLen,:),'color','k','linewidth',3)
%line([modelparams.Tstim_on modelparams.Tstim_on], [0 1.1*max([s.r1; s.r2],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
%line([modelparams.Tstim_off modelparams.Tstim_off], [0 1.1*max([s.r1; s.r2],[],'all')], 'linewidth', 1.5, 'color', [0.7 0.7 0.7])
ylabel('firing rate (/s)', 'fontsize', 24)
xlabel('time in trial (s)', 'fontsize', 24)
set(gca,'tickDir', 'out')
set(gca,'ytick',0:20:max([s.r1;s.r2]));
set(gcf,'color','w');
set(gca,'xtick',[0:0.5:2])
%% Initialization
modelparams = struct( ...
    'gE', 0.2, ...
    'gI', -0.25, ... % cross-inhibition strength [nA] 
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
Icommon = 0.10;
Idiff = 0.001;
s = simulationDMTwoSti(modelparams,0.05,Icommon+Idiff,Icommon-Idiff,8,false);
%% all tests
tau0test = [0.02 0.01 0.005 0.002];
aTest = [250:10:300];
gItest = [-0.4, -0.30]; %gItest = -0.55:0.1:-0.25;
Icommon = 0.12; %Icommon = 0.05:0.05:0.2;  
%Idiff = -0.0625:0.005:0.0625;
Idiff = -0.0225:0.005:0.0225;
sim = 50;
savepath = '/Volumes/ZS256/MCNproject/simulations/';
%% tau and gI
for T = 1:1

    choices = zeros(length(gItest),length(Icommon),length(Idiff));
    RT = cell(length(gItest),length(Icommon),length(Idiff));
    tests = cell(length(gItest), length(Icommon), length(Idiff), sim);
    errors = -sim*ones(length(gItest),length(Icommon),length(Idiff));
    modelparams.tau0 = tau0test(T);
    %% simulation choice and RTs
    for j = 1:length(gItest)
        modelparams.gI = gItest(j);
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
    
    %save([savepath 'tau0=' num2str(tau0test(T)) '.mat'], 'errors', 'choices', 'gItest', 'Icommon', 'Idiff', 'modelparams', 'RT', 'sim', 'tests')
end
%% simulate without choices
Icommon = 0.12;
Idiff = -0.0225:0.005:0.0225;
sim = 20;
%%
tests = cell(length(gItest), length(Icommon), length(Idiff), sim);

for j = 1:length(gItest)
    modelparams.gI = gItest(j);
    for i = 1:length(Icommon)
        for h = 1:length(Idiff)
            for si = 1:sim
                s = simulationDMTwoSti(modelparams,0.05,Icommon(1)-Idiff(h),Icommon(1)+Idiff(h),8,false);
                tests{j,i,h,si} = s;
            end
        end
    end
end
%%
for j = 1:length(gItest)
    figure; hold on;
    for i = 1:length(Icommon)   
        for h = 1:length(Idiff)
            frTemp = tests{j,i,h};
            for si = 1:sim
                subplot(sim,length(Idiff),length(Idiff)*(si-1)+h); hold on;        
                plot(s.t,frTemp(1:length(s.t),si),'color','r');
                plot(s.t,frTemp(length(s.t)+1:end,si),'color','k')
                if si == 1
                    title(['gI=' num2str(gItest(j)) ' Idiff=' num2str(Idiff(h))])
                end
            end
            
        end
    end
end
%% plotting all simulated traces
for j = 1:length(gItest)
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
%%
% with decision making
I1 = 0.01:0.02:0.3;
I2 = [0.01:0.02:0.3]+0.001;
Idiff = -0.08:0.005:0.08;
Idiff = Idiff + 0.001;
Isum = 0.1:0.1:1;
sim = 20;
choices = zeros(length(Isum),length(Idiff));
RT = cell(length(Isum),length(Idiff));
%%
for i = 1:length(Isum)
    for j = 1:length(Idiff)
        for t = 1:sim
            s = simulationDM(modelparams,Isum(i)+Idiff(j),Isum(i)-Idiff(j),8,true);
            choices(i,j) = choices(i,j) + s.choice;
            RT{i,j} = [RT{i,j} s.RT];
        end
    end
end
%% plotting all choices and RTs
for h = 1:length(gItest)
figure;
imagesc(squeeze(choices(h,:,:))/sim);
title(['gI = ' num2str(gItest(h))])
figure;
subplot(1,2,1);
imagesc(cellfun(@mean,squeeze(RT(h,:,:))));
title(['gI = ' num2str(gItest(h))])
subplot(1,2,2);
imagesc(cellfun(@std,squeeze(RT(h,:,:))));
title(['gI = ' num2str(gItest(h))])
end
%%
%% plotting all choices and RTs
for h = 1:2
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

for h = 1:length(gItest)
    figure; hold on;
    colors = cool(size(choices,2));
    currChoices = squeeze(choices(h,:,:));
    for i = 1:size(choices,2)
        p = plot(Idiff, 0.5*(currChoices(i,:)/sim + 1), 'linewidth', 1.5, 'color', colors(i,:));
    end
    legend(cellfun(@num2str,num2cell(Icommon), 'UniformOutput', false))
    title(['p(choice == 1)' 'gI = ' num2str(gItest(h))]);
end
%% plot action selection sep by Idiff & gI
colors = cool(size(choices,1));
figure; hold on;
for h = 1:length(gItest)
    currChoices = squeeze(choices(h,:,:));
    currChoices = mean(currChoices,1);
    plot(2*Idiff, 0.5*(currChoices/sim + 1), 'linewidth', 2.5, 'color', colors(h,:));
end
legend(cellfun(@num2str,num2cell(gItest), 'UniformOutput', false), 'FontSize', 18)
xlabel('I_2 - I_1', 'FontSize', 24)
ylabel('P(choice=2)','FontSize', 24)
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'XTick',[])
set(gca,'ytick',[0 1],'FontSize', 18)
xlim([-0.08, 0.08])
%% scatter action selection sep by Idiff & gI
colors = cool(size(choices,1));
figure; hold on;
for h = 1:length(gItest)
    currChoices = squeeze(choices(h,:,:));
    currChoices = mean(currChoices,1);
    scatter(Idiff, 0.5*(currChoices/sim + 1), 7, colors(h,:));
end
legend(cellfun(@num2str,num2cell(gItest), 'UniformOutput', false))
%% RT with Idiff Icommon







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


%% fit all choice curves
savepath = '/Volumes/ZS256/MCNproject/simulations/';
tau0Test = [0.02, 0.005 0.002];
sigma = zeros(length(tau0Test), length(gItest),3);
for i = 1:length(tau0Test)
    load([savepath 'tau0=' num2str(tau0Test(i)) '.mat'])
    x = Idiff*2;
     for j = 1:length(gItest)
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
    errorbar(gItest,squeeze(sigma(i,:,1)),0.5*(squeeze(sigma(i,:,3)-sigma(i,:,2))), 'LineWidth', 3, 'Color', colors(i,:))
end
legend(cellfun(@num2str,num2cell(tau0Test), 'UniformOutput', false),'FontSize', 18);
xlabel('gI', 'FontSize', 24)
ylabel('\sigma','FontSize', 24)
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'XTick', -0.5:0.1:-0.2,'FontSize', 18)
set(gca,'ytick',0:0.01:0.04,'FontSize', 18)
xlim(minmax(gItest))
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
%% simulate exploration
x = -0.1:0.001:0.1;
y = 1./(1+exp(-x/0.01));
figure; hold on;
plot(x,y,'linewidth', 2.5, 'color', 'k')
% plot(x, 0.5*sign(x)+0.5,'linewidth', 2.5, 'color', 'm')
xlabel('c_2 - c_1', 'FontSize', 24)
ylabel('P(choice=2)','FontSize', 24)
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'XTick',[])
set(gca,'ytick',[0 1],'FontSize', 18)
xlim([-0.08, 0.08])
%% RT plot
RTrecons = RT(:,:,1:(0.5*length(Idiff)));
for i = 1:length(gItest)
    for j = 1:length(Icommon)
        for h = 1:0.5*length(Idiff)
            RTrecons{i,j,h} = [RTrecons{i,j,h}, RT{i,j,length(Idiff)+1-i}];
        end 
    end
end
%%
IdiffSize = flip(sort(unique(abs(Idiff))));
colors = cool(length(Icommon));
figure;
for i = 1:length(gItest)
    subplot(1,length(gItest),i); hold on;
    for j = 1:length(Icommon)
        m = squeeze(cellfun(@mean, RTrecons(i,j,:)));
        err = squeeze(cell2mat(cellfun(@std, RTrecons(i,j,:),'UniformOutput', false)))/sqrt(sim*2-1);
        errorbar(IdiffSize,m,err,'color',colors(j,:),'linewidth',2);
        ylim([0 0.55])
    end
    legend(cellfun(@num2str,num2cell(Icommon), 'UniformOutput', false))
    title(['gI=' num2str(gItest(i))])
end
%%
RTgI = cell(length(gItest),0.5*length(Idiff));
for i = 1:length(gItest)
    for j = 1:0.5*length(Idiff)
        for h = 1:length(Icommon)
            RTgI{i,j} = [RTgI{i,j} RTrecons{i,h,j}];
        end
    end
end
%%
figure;
colors = cool(length(gItest));
hold on;
for i = 2:length(gItest)
    m = squeeze(cellfun(@mean, RTgI(i,:)));
    err = squeeze(cell2mat(cellfun(@std, RTgI(i,:),'UniformOutput', false)))/sqrt(sim*8-1);
    errorbar(IdiffSize,m,err,'color',colors(i,:),'linewidth',2);
end
legend(cellfun(@num2str,num2cell(gItest(2:end)), 'UniformOutput', false))
%%