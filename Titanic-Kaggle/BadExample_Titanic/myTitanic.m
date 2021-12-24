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

writetable(normalized_data,'process.csv')