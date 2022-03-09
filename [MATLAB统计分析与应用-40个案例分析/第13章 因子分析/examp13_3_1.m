%--------------------------------------------------------------------------
%  examp13.3-1  从相关系数矩阵出发进行因子分析
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%***************************定义相关系数矩阵PHO*****************************
PHO = [1     0.79    0.36    0.76    0.25    0.51
          0.79  1       0.31    0.55    0.17    0.35
          0.36  0.31    1       0.35    0.64    0.58
          0.76  0.55    0.35    1       0.16    0.38
          0.25  0.17    0.64    0.16    1       0.63
          0.51  0.35    0.58    0.38    0.63    1
          ];

%******************调用factoran函数根据相关系数矩阵作因子分析*****************
% 从相关系数矩阵出发，进行因子分析（不进行因子旋转）
[lambda,psi,T] = factoran(PHO,2,'xtype','covariance','delta',0,'rotate','none')


% 定义元胞数组，以元胞数组形式显示结果
head = {'变量', '因子f1', '因子f2'};
varname = {'身高','坐高','胸围','手臂长','肋围','腰围','<贡献率>','<累积贡献率>'}';
Contribut = 100*sum(lambda.^2)/6;
CumCont = cumsum(Contribut);
result1 = num2cell([lambda; Contribut; CumCont]);
result1 = [head; varname, result1]


%从相关系数矩阵出发，进行因子分析（进行因子旋转）
[lambda,psi,T] =factoran(PHO,2,'xtype','covariance','delta',0)
Contribut = 100*sum(lambda.^2)/6
CumCont = cumsum(Contribut)


% 从相关系数矩阵出发，进行因子分析，公共因子数为3（进行因子旋转）
[lambda,psi,T] = factoran(PHO,3,'xtype','covariance','delta',0)
Contribut = 100*sum(lambda.^2)/6
CumCont = cumsum(Contribut)