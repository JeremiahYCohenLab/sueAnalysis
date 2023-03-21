function s = behAnalysisNoPlotSim_opMD(allChoices, allRewards, laser)
    allNoRewards = allChoices;
    allNoRewards(allRewards~=0) = 0;
    s.allNoRewards = allNoRewards;
    s.laser = laser;
    s.allChoices = allChoices;
    s.allRewards = allRewards;
    s.rwd_Inds = find(abs(allRewards)>0);
    s.nrwd_Inds = find(abs(allRewards)==0);
    s.changeChoice_Inds = find(allChoices(2:end)~=allChoices(1:end-1))+1;
    s.stayChoice_Inds = find(allChoices(2:end)==allChoices(1:end-1))+1;
end