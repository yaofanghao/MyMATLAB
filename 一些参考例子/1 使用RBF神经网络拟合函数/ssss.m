
%% ѧϰĿ�꣺ ��RBF������ʵ�ַ����Եĺ����ع�
clc 
clear all

%% �����������x1 x2
x1=-1:0.01:1;
x2=-1:0.01:1;

%%  �����������y
y=60+x1.^2-6*cos(6*pi*x1)+6*x2.^2-6*cos(6*pi*x2); 
 
%%  ����RBF����
net=newrbe([x1;x2],y)
 
%%  �������
t=sim(net,[x1;x2]);
 
%%  �������Ч��ͼ
figure(1)
plot3(x1,x2,y,'rd');
hold on;
plot3(x1,x2,t,'b-.');
view(100,25)
title(' RBF����������Ч��')
xlabel('x1')
ylabel('x2')
zlabel('y')
grid on 


%%   ����QQ��1960009019
%%   ���߽���΢�Ź��ںţ�����һƷ�� 
