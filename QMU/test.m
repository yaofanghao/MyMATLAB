% doc random

% fitdist 
% 拟合概率分布

%%  蒙特卡洛法小练习
function piva = PiMonteCarlo(n)
% 用蒙特卡洛模拟法求圆周率pi
% n为投点个数，可以是标量或向量
x = 0;y = 0;d = 0;
m = length(n);
pivalue = zeros(m,1);
for i = 1:m
    x = 2*rand(n(i),1)-1;
    y = 2*rand(n(i),1)-1;
    d = x.^2+y.^2;
    pivalue(i) = 4*sum(d <= 1)/n(i);
end

if nargout == 0
    if m > 1
        plot(n,pivalue,'.')
        h = refline(0,pi);
        set(h,'linewidth',2,'color','k');
        text(1.05*n(end),pi,'\pi','fontsize',15);
        xlabel('投点个数');
        ylabel('\pi的模拟值');
    else
        plot(x,y,'.')
        hold on
        h = rectangle('Position',[-1 -1 2 2],'LineWidth',2);
        t = linspace(0,2*pi,100);
        plot(cos(t),sin(t),'r','linewidth',2);
        xlabel('X');
        ylabel('Y');
        axis([-1.1 1.1 -1.1 1.1]);
        axis equal;
    end
else
    piva = pivalue;
end

%% 从高斯分布N(5, 3)中抽样
% x = norminv(linspace(0+eps, 1, 10000), 5, 3);
x = logninv(linspace(0+eps, 1, 10000), 5, 3);
idx = randperm(length(x));
y = x(idx);
 
figure
plot(y, '.')
 
figure
histogram(y, 'binwidth', 0.1)

%% 
