addpath(genpath(pwd));
warning off;
%% Classification
load C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabLibraries\ROC_Plot\ROC_Plot\Data\Train_dataset.mat;
load C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabLibraries\ROC_Plot\ROC_Plot\Data\Train_label.mat;
load C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabLibraries\ROC_Plot\ROC_Plot\Data\Test_dataset.mat;
load C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabLibraries\ROC_Plot\ROC_Plot\Data\Test_label.mat;
load C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabLibraries\ROC_Plot\ROC_Plot\test_features.mat;

%% K - Nearest neighbor Classifier
[ KNN_Results ] = KNN_Classification( Feature_Vector',Train_Feature_Vector,Train_label,Test_Feature_Vector,Test_label );
