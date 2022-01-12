clc
clear all
tic

%% __________________________________________________________________________________________
%1.抛硬币实验
%0为反面，1为正面，求出现正面的概率
num=[10 100 200 1000 10000  1000000];
for i=1:length(num)
coin=randi(2,num(i),1)-1;
numz=find(coin==1);
pz(i)=length(numz)/length(coin);
end
fprintf('1.分别抛硬币%3.0f  %3.0f %3.0f  %4.0f  %5.0f  %8.0f次\n',num)
fprintf('正面的概率为%3.2f %3.2f %3.2f %4.2f %5.2f %8.2f\n\n',pz)
%% __________________________________________________________________________________________
%2.模拟计算圆周率
x = -1 + 2*rand(1e8,1); %生成【-1 1】之间的随机数
y = -1 + 2*rand(1e8,1);
% 搜寻落到半径为1的圆内的个数  x.^2 + y.^2 <= 1
z = x(x.^2 + y.^2 <= 1);
% 精确解. 3.1416
pitest = 4*length(z)/length(x);  %原理：所有随机数落到半径为1的圆内的概率等于圆的面积比矩形的面积 pi*1^2/4=length(z)/length(x)
err = abs(pitest - pi)/pi;
fprintf('2.梦塔卡罗模拟π的值为%3.4f,相对误差为%3.4f\n\n ',pitest,err)
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
%% 3.要弹簧秤在90%的可靠度区间内工作，弹簧秤最大伸长量为多少？
nn=10000;
m=10+normrnd(0,1,1,nn);%正态分布m随机抽样
k=1000+normrnd(0,50,1,nn);%正态分布k随机抽样
dl=(9.8*m./k-10*9.8/1000)*1000; %d单位mm
[f,xi]=ksdensity(dl);
figure(2)
subplot(1,2,1)
plot(xi,f,'r')
xlabel('弹簧伸长超过正常量')
ylabel('对应超过伸长量出现概率')
title('弹簧伸长超过正常量概率密度分布图')
subplot(1,2,2)
h1 = histogram(dl);
xlabel('弹簧伸长超过正常量')
ylabel('对应弹簧伸长超过正常量出现的次数')
title('弹簧伸长超过正常量直方图')
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
maxl=10*9.8+maxdl; %正常伸长量加上超过正常的伸长量
fprintf('3.可靠度为%3.3f时，弹簧最大伸长量为%3.2f\n',pr,maxl)
fprintf('可靠度误差%3.4f\n',epr)
toc
%__________________________________________________________________________________________