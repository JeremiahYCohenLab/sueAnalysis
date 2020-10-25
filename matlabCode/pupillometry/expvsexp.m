%% fit control and inhi
[root,sep] = currComputer_df;
sessionFolder = 'session';
ani = 'combine';
g1 = 'control';
g2 = 'Inh';
workbookFile = 'Z:\combineAnimals'; 
[A, dayList, C] = xlsread(workbookFile, ani);


[~,col] = find(~cellfun(@isempty,strfind(dayList, g1)) == 1);
dayListG1 = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayListG1),1);
if ~isempty(endInd)
    dayListG1 = dayListG1(1:endInd-1,:);
end

for i = 1:length(dayListG1)
    sprintf(dayListG1{i})
   [~, states] = fitHmm(dayListG1{i},1);
   pEeG1(i) = length(find(states == 3))/length(states);
   
end


[~,col] = find(~cellfun(@isempty,strfind(dayList, g2)) == 1);
dayListG2 = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayListG2),1);
if ~isempty(endInd)
    dayListG2 = dayListG2(1:endInd-1,:);
end

for i = 1:length(dayListG2)
    sprintf(dayListG2{i})
   [~, states] = fitHmm(dayListG2{i},1);
   pEeG2(i) = length(find(states == 3))/length(states);
end
%% fit all sessions
[root,sep] = currComputer;
sessionFolder = 'session';
ani = 'combine';
g = 'Fair';
workbookFile = 'Z:\combineAnimals'; 
[A, dayList, C] = xlsread(workbookFile, ani);


[~,col] = find(~cellfun(@isempty,strfind(dayList, g)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

for i = 1:length(dayList)
%    sprintf(dayList{i})
   [~, states] = fitHmmOpt(dayList{i},1);
end

