animals = {'ZS066','ZS068','ZS069','ZS070','ZS071'};
dates = {'20211027', '20211028', '20211029','20211030', '20211031'};
sessions = cell(length(animals), length(dates));
% generate all session names
for i = 2:length(animals)
    for j = 1:length(dates)
        sessionName = ['m' animals{i} 'd' dates{j}];
%         timeAlign(sessionName, 1);
        pupilSessionAnalysis(sessionName);
    end
end
%%