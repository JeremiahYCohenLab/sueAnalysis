figure2;
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
csOn = [behSessionData(s.CSplus).CSon];
subplot(2, 1, 1); plot(ledLL);
subplot(2, 1, 2); scatter([csOn-csOn(1)]/1000 * 25.5 + 107, ones(size(csOn)));
%%
ledLL(1) = -1;
ledLL(end) = -1;
upEdge = find(ledLL(2:end)>0 & ledLL(1:end-1)<0);
downEdge = find(ledLL(2:end)<0 & ledLL(1:end-1)>0);
