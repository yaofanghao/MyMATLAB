%%  ѧϰĿ�꣺ʹ���Ŵ��㷨��ⲻ��ʽ
clear all
clc
A = [2 1; -1 2; 2 2];
b = [3; 2; 3];
lb = zeros(2,1);
[x,fval,exitflag] = ga(@lincontest6,2,A,b,[],[],lb)
%%   ����QQ��1960009019
%%   ���߽���΢�Ź��ںţ�����һƷ�� 
