clc
clear all
tic

%% __________________________________________________________________________________________
%1.��Ӳ��ʵ��
%0Ϊ���棬1Ϊ���棬���������ĸ���
num=[10 100 200 1000 10000  1000000];
for i=1:length(num)
coin=randi(2,num(i),1)-1;
numz=find(coin==1);
pz(i)=length(numz)/length(coin);
end
fprintf('1.�ֱ���Ӳ��%3.0f  %3.0f %3.0f  %4.0f  %5.0f  %8.0f��\n',num)
fprintf('����ĸ���Ϊ%3.2f %3.2f %3.2f %4.2f %5.2f %8.2f\n\n',pz)
%% __________________________________________________________________________________________
%2.ģ�����Բ����
x = -1 + 2*rand(1e8,1); %���ɡ�-1 1��֮��������
y = -1 + 2*rand(1e8,1);
% ��Ѱ�䵽�뾶Ϊ1��Բ�ڵĸ���  x.^2 + y.^2 <= 1
z = x(x.^2 + y.^2 <= 1);
% ��ȷ��. 3.1416
pitest = 4*length(z)/length(x);  %ԭ������������䵽�뾶Ϊ1��Բ�ڵĸ��ʵ���Բ������Ⱦ��ε���� pi*1^2/4=length(z)/length(x)
err = abs(pitest - pi)/pi;
fprintf('2.��������ģ��е�ֵΪ%3.4f,������Ϊ%3.4f\n\n ',pitest,err)
figure(1)
xlabel('x')
ylabel('y')

hold on
  axis equal
        set(gca,'XLim',[-1 1]);
        set(gca,'YLim',[-1 1]);
for k = 1:100 %length(x)
    if x(k)^2 + y(k)^2 <= 1
        % the point is inside the sphere, show it in red
        plot(x(k), y(k), '.r')
    else
        % the point is outside the sphere, show it in blue
        plot(x(k), y(k) ,'.b') 
    end 
    title('%3.4f',k)
   drawnow update
   %text(0,0,num2str(numindex(1)),'HorizontalAlignment','center','VerticalAlignment','middle');
end
%__________________________________________________________________________________________
%% 3.Ҫ���ɳ���90%�Ŀɿ��������ڹ��������ɳ�����쳤��Ϊ���٣�
nn=10000;
m=10+normrnd(0,1,1,nn);%��̬�ֲ�m�������
k=1000+normrnd(0,50,1,nn);%��̬�ֲ�k�������
dl=(9.8*m./k-10*9.8/1000)*1000; %d��λmm
[f,xi]=ksdensity(dl);
figure(2)
subplot(1,2,1)
plot(xi,f,'r')
xlabel('�����쳤����������')
ylabel('��Ӧ�����쳤�����ָ���')
title('�����쳤���������������ܶȷֲ�ͼ')
subplot(1,2,2)
h1 = histogram(dl);
xlabel('�����쳤����������')
ylabel('��Ӧ�����쳤�������������ֵĴ���')
title('�����쳤����������ֱ��ͼ')
[mua,stda] = normfit(dl);
rfg=norminv(0.9,mua,stda);
maxdl=abs(rfg);
s=0;
for i=1:nn
    if  dl(i)<maxdl 
        s=s+1;
    end
end
pr=s/nn;
epr=abs(0.9-pr)/0.9;
maxl=10*9.8+maxdl; %�����쳤�����ϳ����������쳤��
fprintf('3.�ɿ���Ϊ%3.3fʱ����������쳤��Ϊ%3.2f\n',pr,maxl)
fprintf('�ɿ������%3.4f\n',epr)
toc
%__________________________________________________________________________________________