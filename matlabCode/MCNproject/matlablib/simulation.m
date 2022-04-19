function s = simulation(p,coh,n_trial)
NT = round(p.Ttotal/p.dt);
t_plot = [1:NT] * p.dt;
t_stim = (t_plot > p.Tstim_on) & (t_plot < p.Tstim_off);

mean_stim = ones(1,NT) * p.mu0 * p.Jext/1000; %[nA]
diff_stim = p.Jext * p.mu0 * coh/100. * 2;

Istim1_plot = (mean_stim + diff_stim/2/1000) .* t_stim; %[nA]
Istim2_plot = (mean_stim - diff_stim/2/1000) .* t_stim;

% Initialize S1 and S2
S1 = 0.1 * ones(n_trial,1);
S2 = 0.1 * ones(n_trial,1);

Ieta1 = zeros(n_trial,1);
Ieta2 = zeros(n_trial,1);

n_record = round(p.record_dt/p.dt);
i_record = 1;
N_record = round(p.Ttotal/p.record_dt);
s.r1 = zeros(N_record, n_trial);
s.r2 = zeros(N_record, n_trial);
s.t  = zeros(N_record, 1);
s.I1 = zeros(N_record, 1);
s.I2 = zeros(N_record, 1);
s.S1 = zeros(N_record, n_trial);
s.S2 = zeros(N_record, n_trial);
s.Ieta1 = zeros(N_record, n_trial);
s.Ieta2 = zeros(N_record, n_trial);
% Loop over time points in a trial
for i_t = 1:NT
    % Random dot stimulus
    Istim1 = Istim1_plot(i_t);
    Istim2 = Istim2_plot(i_t);

    % Total synaptic input
    Isyn1 = p.gE * S1 + p.gI * S2 + Istim1 + Ieta1;
    Isyn2 = p.gE * S2 + p.gI * S1 + Istim2 + Ieta2;

    % Transfer function to get firing rate

    r1  = FIcurve(Isyn1);
    r2  = FIcurve(Isyn2);

    %---- Dynamical equations -------------------------------------------

    % Mean NMDA-mediated synaptic dynamics updating
    S1_next = S1 + p.dt * (-S1/p.tauS + (1-S1) .* p.gamma .*r1);
    S2_next = S2 + p.dt * (-S2/p.tauS + (1-S2) .* p.gamma .*r2);

    % Ornstein-Uhlenbeck generation of noise in pop1 and 2
    Ieta1_next = Ieta1 + (p.dt/p.tau0) * (p.I0-Ieta1) + sqrt(p.dt/p.tau0) * p.sigma * normrnd(0,1,n_trial,1);
    Ieta2_next = Ieta2 + (p.dt/p.tau0) * (p.I0-Ieta2) + sqrt(p.dt/p.tau0) * p.sigma * normrnd(0,1,n_trial,1);

    if mod(i_t,n_record) == 1
        s.r1(i_record,:) = r1';
        s.r2(i_record,:) = r2';
        s.I1(i_record) = Istim1;
        s.I2(i_record) = Istim2;
        s.t(i_record) = i_t * p.dt;
        s.S1(i_record,:) = S1';
        s.S2(i_record,:) = S2';
        s.Ieta1(i_record,:) = Ieta1';
        s.Ieta2(i_record,:) = Ieta2';
        
        i_record = i_record+1;
    end
    
    S1 = S1_next;
    S2 = S2_next;
    Ieta1 = Ieta1_next;
    Ieta2 = Ieta2_next;
end