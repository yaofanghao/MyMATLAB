clc;clear;
A1 = [0,0,1300];
A2 = [5000,0,1700];
A3 = [0,5000,1700];
A4 = [5000,5000,1300];
Origin = [0,0,0];
Mat_Di = [760,4550,4500,6300];

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
Z_Alfa = inv(G_a' * inv_Big_Fi * G_a) * G_a' * inv_Big_Fi * H;

cov_Z_Alfa = inv(G_a' * inv_Big_Fi * G_a);

B_once = diag([Z_Alfa(1),Z_Alfa(2),Z_Alfa(3),1/2]);

Big_Fi_Sec = 4 * B_once * cov_Z_Alfa * B_once;
inv_Big_Fi_Sec = inv(Big_Fi_Sec);

G_Sec = [1,0,0;
    0,1,0;
    0,0,1;
    1,1,1;];

H_Sec = [Z_Alfa(1)^2, Z_Alfa(2)^2, Z_Alfa(3)^2, Z_Alfa(4)];
H_Sec = H_Sec';

Z_Alfa_Sec = inv(G_Sec' * inv_Big_Fi_Sec * G_Sec) * G_Sec' * inv_Big_Fi_Sec * H_Sec;

Z_Alfa_Sec = sqrt(Z_Alfa_Sec);