%Titanic Machine Learning from Disaster Compteition from Kaggle
clear
titanic_train = readtable('titanic_train.csv'); %加载数据集
titanic_test = readtable('titanic_test.csv');

titanic_train_missing = sum(ismissing(titanic_train));  %处理缺失值
titanic_test_missing = sum(ismissing(titanic_test));

filled_data = titanic_train ;  % New tables for storing filled data
test_filled_data = titanic_test ;

%----------------训练集
% Calculate mean values and cast double value to integer
%------- Age
mean_age = cast(mean(titanic_train.Age, 'omitnan'),'uint8') ; 
filled_age = fillmissing(titanic_train.Age, 'constant', mean_age);
filled_data.Age = filled_age;
%------- Fare
mean_fare = cast(mean(titanic_train.Fare, 'omitnan'),'uint8') ;
filled_fare = fillmissing(titanic_train.Fare, 'constant', mean_fare);
filled_data.Fare = filled_fare;

%---------------- 测试集
%Set test mean age as training mean age
test_mean_age = mean_age ;
test_mean_fare = mean_fare ;
%------- Age
test_filled_age = fillmissing(titanic_test.Age, 'constant', test_mean_age);
test_filled_data.Age = test_filled_age;
%------- Fare
test_filled_fare = fillmissing(titanic_test.Fare, 'constant', test_mean_fare);
test_filled_data.Fare = test_filled_fare;

%-------------------------------处理分类数据
%---------------- 训练集
%Seperate genders into different columns
filled_data = categorical_data_to_dummy_variables(filled_data, filled_data.Sex); 
filled_data.Sex = [];
%---------------- 测试集
test_filled_data = categorical_data_to_dummy_variables(test_filled_data, test_filled_data.Sex); 
test_filled_data.Sex = [];

%-------------------------------处理异常值(outliers)
plot(filled_data.Age) %Age varies between 0-80 which can be accepted as normal
%Remove rows which has age less than 1
toDelete = filled_data.Age < 1;
filled_data(toDelete,:) = [];
%Remove rows which has age not integer
toDelete2 = mod(filled_data.Age,1) ~= 0;
filled_data(toDelete2,:) = [];

%-------------------------------特征缩放 (归一化)
%---------------- Training set
%New table for normalized data
normalized_data = filled_data ;
%Feature scaling for the Age
normalized_age = (filled_data.Age - min(filled_data.Age)) / (max(filled_data.Age) - min(filled_data.Age));
normalized_data.Age = normalized_age;
%Feature scaling for the Fare
normalized_fare = (filled_data.Fare - min(filled_data.Fare)) / (max(filled_data.Fare) - min(filled_data.Fare));
normalized_data.Fare= normalized_fare;
%Feature scaling for the SibSp
normalized_sibsp = (filled_data.SibSp - min(filled_data.SibSp)) / (max(filled_data.SibSp) - min(filled_data.SibSp));
normalized_data.SibSp = normalized_sibsp;
%Feature scaling for the Parch
normalized_parch = (filled_data.Parch - min(filled_data.Parch)) / (max(filled_data.Parch) - min(filled_data.Parch));
normalized_data.Parch = normalized_parch;
%Feature scaling fot the Pclass
normalized_pclass = (filled_data.Pclass - min(filled_data.Pclass)) / (max(filled_data.Pclass) - min(filled_data.Pclass));
normalized_data.Pclass = normalized_pclass;

%---------------- Test set
%New table for normalized data
test_normalized_data = test_filled_data ;
%Feature scaling for the Age
test_normalized_age = (test_filled_data.Age - min(test_filled_data.Age)) / (max(test_filled_data.Age) - min(test_filled_data.Age));
test_normalized_data.Age = test_normalized_age;
%Feature scaling for the Fare
test_normalized_fare = (test_filled_data.Fare - min(test_filled_data.Fare)) / (max(test_filled_data.Fare) - min(test_filled_data.Fare));
test_normalized_data.Fare= test_normalized_fare;
%Feature scaling for the SibSp
test_normalized_sibsp = (test_filled_data.SibSp - min(test_filled_data.SibSp)) / (max(test_filled_data.SibSp) - min(test_filled_data.SibSp));
test_normalized_data.SibSp = test_normalized_sibsp;
%Feature scaling for the Parch
test_normalized_parch = (test_filled_data.Parch - min(test_filled_data.Parch)) / (max(test_filled_data.Parch) - min(test_filled_data.Parch));
test_normalized_data.Parch = test_normalized_parch;
%Feature scaling fot the Pclass
test_normalized_pclass = (test_filled_data.Pclass - min(test_filled_data.Pclass)) / (max(test_filled_data.Pclass) - min(test_filled_data.Pclass));
test_normalized_data.Pclass = test_normalized_pclass;

%-------------------------------- 分类 (决策树模型)
classification_model = fitctree(normalized_data, 'Survived~Age+Fare+Parch+SibSp+female+male+Pclass'); %Classification Model

%-------------------------------- LOOP FOR GENERAL/AVERAGE ACCURACY
general_accuracy = 0;
for a = 1:1
   %分割训练集
   cv = cvpartition(classification_model.NumObservations,'HoldOut', 0.03); %Built-in function for partitioning
   cross_validated_model = crossval(classification_model, 'cvpartition', cv); %Use training part of training set only to built model 
    
   %预测
   Predictions = predict(cross_validated_model.Trained{1}, normalized_data(test(cv),1:end-1));
   
   %--------------------------------ANALYZING THE RESULT
   %Confusion Matrix: / diagonal will give the false predictions, \ will be the rigth predictions.
   Results = confusionmat(cross_validated_model.Y(test(cv)),Predictions); 
   % coloredResultsMatrix = confusionchart(Results,'DiagonalColor','green');
   right_results = Results(1,1) + Results(2,2);
   wrong_results = Results(1,2) + Results(2,1);
   truth_score = right_results /(right_results + wrong_results);
   %Sum each accuracy to calculate general_accuracy outside the loop
   general_accuracy = general_accuracy + truth_score;
end

%Calculate the general accuracy for a number of iterations
general_accuracy = general_accuracy / a; 
    
%Print general accuracy 
disp('General accuracy is:');
disp(general_accuracy);

%-------------------------------- PERFORM PREDICTIONS ON TEST SET FOR KAGGLE
%Use created classification_model's training set and perform predictions on test_normalized data.
test_predictions = predict(cross_validated_model.Trained{1}, test_normalized_data(1:end,1:end-1));
%---------------- Prepare Table for Kaggle
Resulttable = table(test_normalized_data.PassengerId ,int16(test_predictions));
%Set headers for columns
Resulttable.Properties.VariableNames{1} = 'PassengerId';
Resulttable.Properties.VariableNames{2} = 'Survived';
%Write table into file system
writetable(Resulttable, 'resultedtable.csv');

%-------------------------------- VISUALIZE THE RESULTS FOR TRAINING SET
view(cross_validated_model.Trained{1}, 'Mode', 'Graph');

%-------------------------------- FUNCTION TO HANDLE UNORDERED CATEGORICAL DATA
function data = categorical_data_to_dummy_variables(data,variable)
unique_values = unique(variable);
for i=1:length(unique_values)
    dummy_variable(:,i) = double(ismember(variable,unique_values{i})) ;
end 
T = table;
[rows, col] = size(dummy_variable);
for i=1:col
    T1 = table(dummy_variable(:,i));
    T1.Properties.VariableNames = unique_values(i);
    T = [T T1];
end 
    data = [T data]; 
    end