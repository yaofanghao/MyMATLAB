%% ѧϰĿ�꣺ʹ�ý����õ������磨ѵ���ò����棬�´�ֱ�ӵ��ø������磩���з���

clear all;
close all;
P=[-0.4 -0.4 0.5 -0.2 -0.7;-0.6 0.6 -0.4 0.3 0.8];      %��������
T=[1 1 0 0 1];                                          %�������
plotpv(P,T);                                            %��������
net=newp(minmax(P),1,'hardlim','learnpn');              %����������
hold on;
linehandle=plot(net.IW{1},net.b{1});
E=1;
net.adaptParam.passes=10;
while mae(E)                                            %���ﵽҪ���ֹͣѵ��
    [net,Y,E]=adapt(net,P,T);                           %���и�֪���������ѵ��
    linehandle=plotpc(net.IW{1},net.b{1},linehandle);
    drawnow;
end
save net1 net;                                          %��ѵ���õ���������б���
set(gcf,'position',[60,60,300,300]);

%%  �øղŽ�������������з���
clear all;
close all;
load net1.mat;                                  %�����ϴ�ѵ���õ�������
X=[-0.3 0.3 0.9;-0.6 0.2 0.8];                  %��������
Y=sim(net,X);                                   %��������з���
figure;
plotpv(X,Y);                                    %����������
plotpc(net.IW{1},net.b{1});                     %���Ʒ�����
set(gcf,'position',[60,60,300,300]);

%%  