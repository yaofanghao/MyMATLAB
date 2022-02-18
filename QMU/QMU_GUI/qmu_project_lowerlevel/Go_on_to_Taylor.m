%check Res
function [X,Res_All] = Go_on_to_Taylor(d)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明

A1 = [0,0,1300];
A2 = [5000,0,1700];
A3 = [0,5000,1700];
A4 = [5000,5000,1300];
Origin = [0,0,0];
Mat_Di = d;

K(1) = Cal_Dis_3Dim(A1,Origin)^2;
K(2) = Cal_Dis_3Dim(A2,Origin)^2;
K(3) = Cal_Dis_3Dim(A3,Origin)^2;
K(4) = Cal_Dis_3Dim(A4,Origin)^2;

ri_Matrix = Mat_Di.^2;

H = ri_Matrix - K;
H = H';

G_a = [-2 *A1(1), -2 *A1(2), -2 *A1(3), 1;
    -2 *A2(1), -2 *A2(2), -2 *A2(3), 1;
    -2 *A3(1), -2 *A3(2), -2 *A3(3), 1;
    -2 *A4(1), -2 *A4(2), -2 *A4(3), 1;];

B = diag(Mat_Di);
Q = eye(4);

Big_Fi = 4 * B * Q * B;

inv_Big_Fi = inv(Big_Fi);
% 计算 Z_Alfa_1 用于下一步求Z_Alfa_2
Z_Alfa_1 = inv(G_a' * G_a) * G_a'  * H;

% Z_Alfa_1 中三个坐标值
X = sqrt(Z_Alfa_1(1:3));

% 计算跟锚点距离 计算D_i0

Di0_mao(1) = Cal_Dis_3Dim(X,A1);
Di0_mao(2) = Cal_Dis_3Dim(X,A2);
Di0_mao(3) = Cal_Dis_3Dim(X,A3);
Di0_mao(4) = Cal_Dis_3Dim(X,A4);

% 用于计算协方差矩阵Q_Improved
Sigma_Matirx = abs(d - Di0_mao);

B_Improved = diag(Di0_mao);

Q_Improved = diag(Sigma_Matirx).^2;

Big_Fi_Improved = 4 * B_Improved * Q_Improved * B_Improved;

inv_Big_Fi_Improved = inv(Big_Fi_Improved);
Mid_Cal = inv(G_a' * inv_Big_Fi_Improved * G_a);
Z_Alfa_2 = Mid_Cal * G_a' * inv_Big_Fi_Improved * H;

X_Pos_2 = Z_Alfa_2(1:3);
X_Pos_2 = abs(X_Pos_2);
Res(1) = Cal_Dis_3Dim(X_Pos_2,A1);
Res(2) = Cal_Dis_3Dim(X_Pos_2,A2);
Res(3) = Cal_Dis_3Dim(X_Pos_2,A3);
Res(4) = Cal_Dis_3Dim(X_Pos_2,A4);

% 计算Res平方和 （RSS）
Res_All = (Res - Mat_Di).^2;
Res_All = sum(Res_All);


%% Improved Algo
%Z_Alfa = inv(G_a' * inv_Big_Fi * G_a) * G_a' * inv_Big_Fi * H;

% inv_Big_Fi = inv(Big_Fi);


%%
% Z_Alfa = inv(G_a' * inv_Big_Fi * G_a) * G_a' * inv_Big_Fi * H;
% cov_Z_Alfa = inv(G_a' * inv_Big_Fi * G_a);
% 
% B_once = diag([Z_Alfa(1),Z_Alfa(2),Z_Alfa(3),1/2]);
% 
% Big_Fi_Sec = 4 * B_once * cov_Z_Alfa * B_once;
% inv_Big_Fi_Sec = inv(Big_Fi_Sec);
% 
% G_Sec = [1,0,0;
%     0,1,0;
%     0,0,1;
%     1,1,1;];
% 
% H_Sec = [Z_Alfa(1)^2, Z_Alfa(2)^2, Z_Alfa(3)^2, Z_Alfa(4)];
% H_Sec = H_Sec';
% 
% Z_Alfa_Sec = inv(G_Sec' * inv_Big_Fi_Sec * G_Sec) * G_Sec' * inv_Big_Fi_Sec * H_Sec;
% 
% Z_Alfa_Sec = sqrt(abs(Z_Alfa_Sec));
% X = Z_Alfa_Sec;
% 
% % Taylor Af Chan
% 
% 
X = X_Pos_2;
cnt = 0;
while 1
    cnt = cnt + 1;
    Di_mao(1) = Cal_Dis_3Dim(X,A1);
    Di_mao(2) = Cal_Dis_3Dim(X,A2);
    Di_mao(3) = Cal_Dis_3Dim(X,A3);
    Di_mao(4) = Cal_Dis_3Dim(X,A4);
    
    Delt_H_t = Mat_Di - Di_mao;
    Delt_H_t = Delt_H_t';
    
    Gt_Taylor = [(X(1)-A1(1))/Di_mao(1), (X(2)-A1(2))/Di_mao(1), (X(3)-A1(3))/Di_mao(1);
        (X(1)-A2(1))/Di_mao(2), (X(2)-A2(2))/Di_mao(2), (X(3)-A2(3))/Di_mao(2);
        (X(1)-A3(1))/Di_mao(3), (X(2)-A3(2))/Di_mao(3), (X(3)-A3(3))/Di_mao(3);
        (X(1)-A4(1))/Di_mao(4), (X(2)-A4(2))/Di_mao(4), (X(3)-A4(3))/Di_mao(4);];
    
    inv_Var_Q_Taylor = inv(Q_Improved);
    
    Min_Cal_Delta = inv(Gt_Taylor' * inv_Var_Q_Taylor * Gt_Taylor);
    
    Dleta_Taylor = Min_Cal_Delta * Gt_Taylor' * inv_Var_Q_Taylor * Delt_H_t;
    
    Delta_Taylor_Dis = sqrt(abs(Dleta_Taylor(1)^2 + Dleta_Taylor(2)^2 + Dleta_Taylor(3)^2));
    
    
    X = X + Dleta_Taylor;
    if Delta_Taylor_Dis < 100 || cnt > 500
        break;
    end
    
end
Res(1) = Cal_Dis_3Dim(X,A1);
Res(2) = Cal_Dis_3Dim(X,A2);
Res(3) = Cal_Dis_3Dim(X,A3);
Res(4) = Cal_Dis_3Dim(X,A4);

% 计算Res平方和 （RSS）
Res_All = (Res - Mat_Di).^2;
Res_All = sum(Res_All);
end

