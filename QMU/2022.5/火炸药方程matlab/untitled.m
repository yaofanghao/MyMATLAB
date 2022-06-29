%% 原始数据 12个
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

%% 第一幅图
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, 44, G1, G2, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, 14, G1, G2, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, 50, G1, G2, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, 100, G1, G2, a,b,c,d,e,g,x,y,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;
xlabel('lambda');ylabel('result')
hold off
legend({'I=44','I=14','I=50','I=100'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第二幅图
I=44;
G1=111;
G2=200;
a=0.01;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,0.222,c,d,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,0,c,d,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,0.667,c,d,e,g,x,y,z);
end

plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'b=0.222','b=0','b=0.667'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第三幅图
I=44;
G1=111;
G2=200;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, 0.01,b,c,d,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, 0,b,c,d,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, 0.23,b,c,d,e,g,x,y,z);
end

plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'a=0.01','a=0','a=0.23'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第四幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,4,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,7,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,15,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,20,y,z);
end
plot(lambda,test1,'r');
hold on;
plot(lambda,test2,'b');
hold on;
plot(lambda,test3,'g');
hold on;
plot(lambda,test4,'*');
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'x=4','x=7','x=15','x=20'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第五幅图
I=44;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, 111, G2, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, 3.1, G2, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, 480, G2, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, 850, G2, a,b,c,d,e,g,x,y,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'G1=111','G1=3.1','G1=480','G1=850'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第六幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,0.222,d,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,0,d,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,0.667,d,e,g,x,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, G2, a,b,1,d,e,g,x,y,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'c=0.222','c=0','c=0.667','c=1'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第七幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,c,0.667,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,c,0,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,c,0.111,e,g,x,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, G2, a,b,c,0.333,e,g,x,y,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'d=0.667','d=0','d=0.111','d=0.333'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第八幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,1.66,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,1,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,2,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,3.2,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'y=1.66','y=1','y=2','y=3.2'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第九幅图
I=44;
G1=111;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, 200, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, 0, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, 100, a,b,c,d,e,g,x,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, 850, a,b,c,d,e,g,x,y,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'G2=200','G2=0','G2=100','G2=850'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第十幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
g=0.667;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);
for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,c,d,0.333,g,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,c,d,0,g,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,c,d,0.222,g,x,y,z);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, G2, a,b,c,d,0.667,g,x,y,z);
end
plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'e=0.333','e=0','e=0.222','e=0.667'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第十一幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
x=4;
y=1.66;
z=2;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);

for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,0.667,x,y,z);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,0,x,y,z);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,1,x,y,z);
end

plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'g=0.667','g=0','g=1'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')

%% 第十二幅图
I=44;
G1=111;
G2=200;
a=0.01;
b=0.222;
c=0.222;
d=0.667;
e=0.333;
g=0.667;
x=4;
y=1.66;

lambda = 0.01:0.01:1;
num=100;
test1=zeros(1,num);
test2=zeros(1,num);
test3=zeros(1,num);
test4=zeros(1,num);

for i=1:100
    test1(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,y,2);
end
for i=1:100
    test2(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,y,1);
end
for i=1:100
    test3(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,y,1.6);
end
for i=1:100
    test4(i)=output(0.01*i, I, G1, G2, a,b,c,d,e,g,x,y,3);
end

plot(lambda,test1);
hold on;
plot(lambda,test2);
hold on;
plot(lambda,test3);
hold on;
plot(lambda,test4);
hold on;

xlabel('lambda');ylabel('result')
hold off
legend({'z=2','z=1','z=1.6','z=3'},'Location','southeast')
xlabel('λ');ylabel('dλ/dt')