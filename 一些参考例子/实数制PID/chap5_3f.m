function [Kpidi,BsJ]=pid_gaf(Kpidi,BsJ)%��������������������������Kpidi,BsJ������ֵ
%����ȫ�ֱ���rin yout timef
global rin yout timef

ts=0.001;
%�������ݺ�����ת��Ϊ���ݺ�����
sys=tf(400,[1,50,0]);
%c2d����:�������������һ����ױ�������������ʱ���״̬�ռ�ģ��ת����ɢʱ��״̬�ռ�ģ�͡�
dsys=c2d(sys,ts,'z');
[num,den]=tfdata(dsys,'v');

rin=1.0;
u_1=0.0;u_2=0.0;
y_1=0.0;y_2=0.0;
x=[0,0,0]';
B=0;
error_1=0;
tu=1;
s=0;
P=100;

for k=1:1:P
   timef(k)=k*ts;
   r(k)=rin;
   %ѭ����ÿһ��u��ֵ
   u(k)=Kpidi(1)*x(1)+Kpidi(2)*x(2)+Kpidi(3)*x(3); 
   %�Ծ���u��ÿһ����ֵ�����޷�
   if u(k)>=10
      u(k)=10;
   end
   if u(k)<=-10
      u(k)=-10;
   end   
   
   yout(k)=-den(2)*y_1-den(3)*y_2+num(2)*u_1+num(3)*u_2;
   error(k)=r(k)-yout(k);
%------------ Return of PID parameters -------------
   u_2=u_1;u_1=u(k);
   y_2=y_1;y_1=yout(k);
   
   x(1)=error(k);                % ���� P
   x(2)=(error(k)-error_1)/ts;   % ���� D
   x(3)=x(3)+error(k)*ts;        % ���� I
   
   error_2=error_1;
   error_1=error(k);
    if s==0
       if yout(k)>0.95 && yout(k)<1.05
          tu=timef(k);
          s=1;
       end 
    end
end

for i=1:1:P
    %abs������������ȡ����ֵ
   Ji(i)=0.999*abs(error(i))+0.01*u(i)^2*0.1;
   B=B+Ji(i);   
   if i>1   
       erry(i)=yout(i)-yout(i-1);
       if erry(i)<0
          B=B+100*abs(erry(i));
       end
   end
end
%������������BsJ��ֵ
BsJ=B+0.2*tu*10;