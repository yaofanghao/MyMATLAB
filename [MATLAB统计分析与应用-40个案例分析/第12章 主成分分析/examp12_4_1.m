%--------------------------------------------------------------------------
%  examp12.4-1  读取examp12_4_1.xls中数据，进行主成分分析
%--------------------------------------------------------------------------

%**************************读取数据，并进行标准化变换************************
[X,textdata] = xlsread('examp12_4_1.xls');
XZ = zscore(X);


%**********************************主成分分析*******************************
% 调用princomp函数根据标准化后原始样本观测数据作主成分分析
[COEFF,SCORE,latent,tsquare] = princomp(XZ)

% 为了直观，定义元胞数组result1，用来存放特征值、贡献率和累积贡献率等数据
explained = 100*latent/sum(latent);
[m, n] = size(X);
result1 = cell(n+1, 4);
result1(1,:) = {'特征值', '差值', '贡献率', '累积贡献率'};
result1(2:end,1) = num2cell(latent);
result1(2:end-1,2) = num2cell(-diff(latent));
result1(2:end,3:4) = num2cell([explained, cumsum(explained)])

% 为了直观，定义元胞数组result2，用来存放前2个主成分表达式的系数数据
varname = textdata(3,2:end)';
result2 = cell(n+1, 3);
result2(1,:) = {'标准化变量', '特征向量t1', '特征向量t2'};
result2(2:end, 1) = varname;
result2(2:end, 2:end) = num2cell(COEFF(:,1:2))

% 为了直观，定义元胞数组result3，用来存放每一个地区总的消费性支出，以及前2个主成分的得分数据
cityname = textdata(4:end,1);
sumXZ = sum(XZ,2);
[s1, id] = sortrows(SCORE,1);
result3 = cell(m+1, 4);
result3(1,:) = {'地区', '总支出', '第一主成分得分y1', '第二主成分得分y2'}; 
result3(2:end, 1) = cityname(id);
result3(2:end, 2:end) = num2cell([sumXZ(id), s1(:,1:2)])

% 为了直观，定义元胞数组result4，用来存放前2个主成分的得分数据，以及(食品+其他)-(衣着+医疗)
cloth = sum(XZ(:,[1,8]),2) - sum(XZ(:,[2,7]),2);
[s2, id] = sortrows(SCORE,2);
result4 = cell(m+1, 4);
result4(1,:) = {'地区','第一主成分得分y1','第二主成分得分y2' ,'(食+其他)-(衣+医)'};
result4(2:end, 1) = cityname(id);
result4(2:end, 2:end) = num2cell([s2(:,1:2), cloth(id)])


%***************************前两个主成分得分散点图***************************
plot(SCORE(:,1),SCORE(:,2),'ko');
xlabel('第一主成分得分');
ylabel('第二主成分得分');
gname(cityname);


%**********************根据霍特林T2统计量寻找极端数据************************
result5 = sortrows([cityname, num2cell(tsquare)],2);
[{'地区', '霍特林T^2统计量'}; result5]


%**************************调用pcares函数重建观测数据************************
for i = 1 : 8
     residuals = pcares(X, i);
     Rate = residuals./X;
     E1(i) = sqrt(mean(residuals(:).^2));
     E2(i) = sqrt(mean(Rate(:).^2));
end
E1
E2