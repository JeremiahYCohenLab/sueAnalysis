session = 'mZS078d20220610';
[root,sep] = currComputer_df;

% convert avi to mat, pupilStruct
pupilStruct = pupil_aviToMat(session);
% get p_mean & i_mean, fit
fits = adjustDR_pupil(pupilStruct);
% Ellipse fit
fits = getPupilDiameter_wrapper(pupilStruct, fits);
% hline
fits = get_hline_pupil(pupilStruct, fits);
% clean fit
fits = cleanPupilFits(pupilStruct, fits);
% time alignment
fits = generatePupilTime_wrapper(session, root, pupilStruct, fits);
 

% Go to C:\Users\zhixiao\Documents\gitRepositories\zhixiaoAnalysis\matlabCode(bab)\currentProjects\pupillometry\scripts




