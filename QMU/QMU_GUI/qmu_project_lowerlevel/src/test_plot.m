Xp_upper = 24.4090;Xp_lower = 21.9372;
Xpu_upper = 25.9012 ;Xpu_lower = 22.9167;
Xth = 24.4;
X0 = 26;
x_miritx = [Xp_upper,Xpu_upper,X0];
% xlim([0 10])
% ylim([-0.4 0.8])
x = 1:10;
y = x_miritx(1)+zeros(10,1);
plot(x,y,'--r'); %蓝色线
hold on;
y = x_miritx(2)+zeros(10,1);
plot(x,y,'-b') ;%蓝色线
hold on;
y = x_miritx(3)+zeros(10,1);
p = plot(x,y,'-b'); %蓝色线
hold on;
rectangle('Position',[1 2 5 6])
rymax = 1.01*max(x_miritx);
rymin = 0.95*min(x_miritx);

axis([0 11 rymin rymax]);
annotation('doublearrow',[0.95/11,0.95/11],[24.45/26.5,25.9/26.5],'LineStyle','-',...
    'color','k','LineWidth',2);