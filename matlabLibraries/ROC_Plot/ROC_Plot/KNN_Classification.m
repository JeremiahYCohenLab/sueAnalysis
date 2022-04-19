function [ KNN_Results] = KNN_Classification( Test_Features,Trainning_Features,Labels,Test_Feature_Vector,Test_label )

Class = fitcknn(Test_Features,Trainning_Features,Labels);
disp('KNN Classification:')
if Class==1
    disp('Identified: Class1');
    result = 'Class1';
else if Class==2
        disp('Identified: Class2 ');
        result = 'Class2';
    else if Class==3
            disp('Identified: Class3 ');
            result = 'Class3';
        end
    end
end
Result = result;

[KNN_Results] = KNN_Classification_Performance( Trainning_Features,Labels,Test_Feature_Vector,Test_label );

disp('Accuracy: ');disp((KNN_Results.accuracy).*100);
disp('Sensitivity: ');disp((KNN_Results.sensitivity).*100);
disp('Specificity: ');disp((KNN_Results.specificity).*100);

KNN_Results.Result = Result;
end


function [KNN_Results] = KNN_Classification_Performance( Trainning_Features,Labels,Test_Feature_Vector,Test_label )
% L1 =length(Trainning_Features(:,1));
% 
% Train = Trainning_Features(1:(L1*0.7),:);
% Test = Trainning_Features(((L1*0.7)+1):end,:);
% Train_Labels = Labels(1:(L1*0.7));
% Test_labels = Labels(((L1*0.7)+1):end);

Train = Trainning_Features;
Test = Test_Feature_Vector;
Train_Labels = Labels;
Test_labels = Test_label;

Predicted_Label = knnclassify(Test,Train,Train_Labels);


stats = confusionmatStats(Test_labels,Predicted_Label);

KNN_Results = stats;

name_class = {'Class1';'Class2';'Class3'};
figure('name','KNN: Confusion Matrix');
draw_cm(stats.confusionMat,name_class,max(unique(Labels)));

cols = [200 45 43; 37 64 180; 0 176 80; 0 0 0]/255;
figure('name','KNN: Precision_Recall Curve');
hold on;
plot(sort(stats.recall,'descend'), sort(stats.precision), '-o', 'color', cols(1,:), 'linewidth', 2)
hold off;
xlabel('Recall'); ylabel('Precision');
end

