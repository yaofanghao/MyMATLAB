clc
clear all
%%
%1.matlab����ͳ�Ʒ������ú���
%��1��.����������ֵ
% mean(X)     % XΪ�������򷵻�X�ľ�ֵ
% mean(A)     % AΪ�����򷵻�ÿ�еľ�ֵ
% mean(A,2)   % AΪ�����򷵻�ÿ�еľ�
% 2������������
% var(X)       % XΪ�������򷵻�X�ķ���
% var(A,[],1)  % AΪ�����򷵻�ÿ�еķ���
% var(A,[],2)  % AΪ�����򷵻�ÿ�еķ���
% 3�����׼��
% std(X)         % XΪ�������򷵻�X�ı�׼��
% std(A,flag,1)  % AΪ�����򷵻�ÿ�еı�׼��(flagȡ0��1��ȡ0��ʾ���ɶ�=��������-1��ȡ1ʱ��ʾ���ɶ�=��������)
% std(A,flag,2)  % AΪ�����򷵻�ÿ�еı�׼��(flagȡ0��1��ȡ0��ʾ���ɶ�=��������-1��ȡ1ʱ��ʾ���ɶ�=��������)
% %% ��4��.��Э����
% cov(X)       % XΪ�������򷵻�X��Э����
% var(A)       % AΪ�����򷵻ظþ����Э���������Խ���Ԫ��Ϊԭ����A���еķ���
% cov(X,Y)     % ���ȳ�������X��Y��Э����
%��5��.�����ϵ��
% corrcoef(X,Y)   % �ȳ�������X��Y�����ϵ��
% corrcoef(A)     % ����A�������������ϵ������
%%
%2.��ͬ����ģ�ͳ�������
% normrnd ��������һ����ֵ�ͱ�׼�����̬�ֲ�
% gamrnd ��������gamma�ֲ���α���������
% chi2rnd �������ɿ����ֲ���α���������
% trnd ��������t�ֲ���α���������
% frnd ��������f�ֲ���α���������
% raylrnd   ��������rayleigh�ֲ���α���������
% lognrnd  �����ֲ��������

nn=5000; %�����������
mu1=10; std1=1;  %������ֵ�ͷ���
dataa=normrnd(mu1,std1,1,nn); %��̬�ֲ�����
[mua,stda] = normfit(dataa)%��̬�ֲ���������
fprintf('���۾�ֵ10,������ֵ%4.2f ���۷���1,��������%4.2f\n', mua,stda)
%%
%3.�������������ܶȺ����ֲ����ۻ������ֲ�
X=[min(dataa):0.1:max(dataa)];
Y = normpdf(X,mu1,std1);%�����ܶȺ���
p = normcdf(X,mu1,std1);%�ۻ��ֲ�����
% figure(1)
subplot(1,2,1)
plot(X,Y,'r')
xlabel('����ֵ')
ylabel('��Ӧ����ֵ���ָ���')
title('�������������ܶȺ����ֲ�')
% figure(2)
subplot(1,2,2)
plot(X,p,'r')
xlabel('����ֵ')
ylabel('��Ӧ����ֵ�ۻ�����')
title('���������ۻ������ܶȷֲ�')
%%
%4.��ɢ��������ĸ��ʷֲ�ͼ��ֱ��ͼ
datab=[1 2 4 2 5 3 8 6 5 7 6 3 4 3 4 9 7  3 4 6 5 5 6 4 3 2 5 7  6 5 5 6 6 4 4 4 5 5 1 3 5 ];
mub=mean(datab);
stdb=std(datab);
fprintf('��ɢ������ֵΪ%4.2f ����Ϊ%4.2f\n', mub,stdb)
[f,xi]=ksdensity(datab);
figure(3)
plot(xi,f,'r')
xlabel('����ֵ')
ylabel('��Ӧ����ֵ���ָ���')
title('������ݸ����ܶȷֲ�ͼ')
figure(4)
h1 = histogram(datab);
xlabel('����ֵ')
ylabel('��Ӧ����ֵ���ֵĴ���')
title('������ݸ���ֱ��ͼ')
%%
%5.���ʷֲ�ģ���������
%dfittool  %���ʷֲ����
save tempresult pd
load  tempresult
rfg=norminv(0.9,pd.mu,pd.sigma);%�ɿ�ָ��
fprintf('����С��%4.1f�ĸ���Ϊ0.9 \n',rfg)
rfd=norminv(0.1,pd.mu,pd.sigma);%�ɿ�ָ��
fprintf('���ݴ���%4.1f�ĸ���Ϊ0.9 \n\n',rfd)
p = normcdf([rfg rfd],pd.mu,pd.sigma)%���ۻ��ֲ�������֤
X1=[min(datab)-2:0.1:max(datab)+2];
Y1 = normpdf(X1,pd.mu,pd.sigma);%�����ܶȺ���
figure(5)
plot(xi,f,'r',X1,Y1,'b.')
xlabel('����ֵ')
ylabel('��Ӧ����ֵ���ָ���')
title('������ݸ����ܶȷֲ�ͼ')
legend('ʵ�ʸ��ʷֲ�','��ϸ��ʷֲ�')
%%