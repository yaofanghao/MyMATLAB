%%   ��ջ��������������
clc;     
clear;
close all;
%%  ����ȫ�ֱ���A
global A
%%   ���복��ͼƬ����ʾ
I=imread('1.png');
figure(1),imshow(I);title('ԭͼ');
%%  
        if length(size(I))==3
            I1 = rgb2gray(I);    %ת��Ϊ�Ҷ�ͼ
        else
           I1 =I;
        end
figure(2),subplot(1,2,1),imshow(I1);title('�Ҷ�ͼ');
figure(2),subplot(1,2,2),imhist(I1);title('�Ҷ�ͼֱ��ͼ');
%%  canny���ӱ�Ե���
I2=edge(I1,'canny',[0.2,0.55]);
figure(3),imshow(I2);title('canny���ӱ�Ե���')
%%  ��ʴ
se=[1;1;1];
I3=imerode(I2,se);         
figure(4),imshow(I3);title('��ʴ��ͼ��');
%%  ƽ��
se=strel('rectangle',[40,40]);
I4=imclose(I3,se);
figure(5),imshow(I4);title('ƽ��ͼ�������');
%% �Ӷ������Ƴ�С���� 
I5=bwareaopen(I4,2000);
figure(6),imshow(I5);
title('�Ӷ������Ƴ�С����');
[y,x,z]=size(I5);
myI=double(I5);
%%   begin����ɨ��
     Blue_y=zeros(y,1);
      for i=1:y
         for j=1:x
             if(myI(i,j,1)==1) 
          %���myI(i,j,1)��myIͼ��������Ϊ(i,j)�ĵ�Ϊ��ɫ
          %��Blue_y����Ӧ�е�Ԫ��white_y(i,1)ֵ��1
           Blue_y(i,1)= Blue_y(i,1)+1;         %��ɫ���ص�ͳ�� 
             end  
         end    
      end
[temp MaxY]=max(Blue_y);
%tempΪ����white_y��Ԫ���е����ֵ��MaxYΪ��ֵ�������� �������е�λ�ã�
PY1=MaxY;
     while ((Blue_y(PY1,1)>=50)&&(PY1>1))
          PY1=PY1-1;
     end    
 PY2=MaxY;
     while ((Blue_y(PY2,1)>=10)&&(PY2<y))
        PY2=PY2+1;
     end
     IY=I(PY1:PY2,:,:);
%IYΪԭʼͼ��I�н�ȡ����������PY1��PY2֮��Ĳ���
 %end����ɨ��
 %%   begin����ɨ��
      Blue_x=zeros(1,x);   %��һ��ȷ��x����ĳ�������
 for j=1:x
         for i=PY1:PY2
             if(myI(i,j,1)==1)
                  Blue_x(1,j)= Blue_x(1,j)+1;               
             end  
         end       
 end
 PX1=1;
      while ((Blue_x(1,PX1)<3)&&(PX1<x))
          PX1=PX1+1;
      end    
 PX2=x;
      while ((Blue_x(1,PX2)<3)&&(PX2>PX1))
             PX2=PX2-1;
      end 
 PX1=PX1-1;
 PX2=PX2+15;
 PY1=PY1-2.5;
 %end����ɨ��
 dw=I(PY1:PY2,PX1:PX2,:);
 imwrite(dw,'dw.jpg'); %��ͼ������д��ͼ����
figure(7),subplot(1,2,1),imshow(IY),title('�з����������');
figure(7),subplot(1,2,2),imshow(dw),title('��λ���к�Ĳ�ɫ����ͼ��')
%%   ����radon�任��ˮƽ����
bw=rgb2gray(dw);
bw1=edge(bw,'sobel','horizontal');
theta=0:179;
r=radon(bw1,theta);
[m,n]=size(r);
c=1;
for i=1:m
    for j=1:n
        if r(1,1)<r(i,j)
            r(1,1)=r(i,j);
            c=j;
        end
    end
end
rot=90-c+2;
pic=imrotate(bw,rot,'crop');
figure(8),subplot(3,1,1),imshow(bw),title('1.��λ������ƻҶ�ͼ��');


figure(8),subplot(3,1,2),imshow(pic),title('2.����radon������ˮƽ�������');
%%  ����radon�任����ֱ����Ľ���
binaryImage = edge(pic,'canny'); 
binaryImage = bwmorph(binaryImage,'thicken'); 
theta = -90:89;
[R,xp] = radon(binaryImage,theta);

[R1,r_max] = max(R);
theta_max = 90;
while (theta_max > 50 || theta_max<-50)
    [R2,theta_max] = max(R1);                      
    R1(theta_max) = 0; 
    theta_max = theta_max - 91;
end
%�Ƕȼ������
H=[1,0,0; tan(-theta_max),1,0;0,0,1];
T=maketform('affine',H);
pic=imtransform(pic,T);
figure(8),subplot(3,1,3), imshow(pic);title('3.����radon��������ֱ�������');


%%   �ַ��ָ�ǰ��Ԥ����
g_max=double(max(max(pic)));
g_min=double(min(min(pic)));
T=round(g_max-(g_max-g_min)/3); % T Ϊ��ֵ������ֵ
[m,n]=size(pic);
d=im2bw(pic,T/256);              % d:��ֵͼ��
figure(9);subplot(2,2,1),imshow(d),title('1.���ƶ�ֵͼ��')
%%   �˲�
h=fspecial('average',3);
d=im2bw(round(filter2(h,d)));
figure(9),subplot(2,2,2),imshow(d),title('2.��ֵ�˲���')
%%   ĳЩͼ����в���
%%   ���ͻ�ʴ
se=strel('square',3); % ʹ��һ��3X3�������ν��Ԫ�ض���Դ�����ͼ������
se=eye(2);        % eye(n) returns the n-by-n identity matrix ��λ����
[m,n]=size(d);
if bwarea(d)/m/n>=0.365
       d=imerode(d,se);
elseif bwarea(d)/m/n<=0.235
       d=imdilate(d,se);
end
figure(9),subplot(2,2,3),imshow(d),title('3.���ͻ�ʴ�����')

%%  ��ͼ��ı߿���вü���ֻ������Ч�ַ�����
[y1,x1,z1]=size(d);
    
I3=double(d);
TT=1;
Y1=zeros(y1,1);
 for i=1:y1
    for j=1:x1
             if(I3(i,j,1)==1) 
                Y1(i,1)= Y1(i,1)+1 ;
            end  
     end       
 end
Py1=1;
Py0=1;
while ((Y1(Py0,1)<30)&&(Py0<y1))
      Py0=Py0+1;
end
Py1=Py0;
 while((Y1(Py1,1)>=30)&&(Py1<y1))
         Py1=Py1+1;
 end
d=d(Py0:Py1,:,:);
figure(9),subplot(2,2,4);
imshow(d),title('4.Ŀ�공������');
%%   ���з����Ͻ��лҶ�ֵ�ۼ�
X1=zeros(1,x1);
for j=1:x1
    for i=1:y1
             if(I3(i,j,1)==1) 
                X1(1,j)= X1(1,j)+1;
            end  
     end       
end
figure(10);
plot(0:x1-1,X1),title('�з������ص�Ҷ�ֵ�ۼƺ�'),xlabel('��ֵ'),ylabel('�ۼ�������');
close all;

Px0=1;
Px1=1;
%%   �ָ��ַ�
for i=1:7
  while ((X1(1,Px0)<2)&&(Px0<x1))
      Px0=Px0+1;
  end
  Px1=Px0;
  while (((X1(1,Px1)>=3)&&(Px1<x1))||((Px1-Px0)<10))
      Px1=Px1+1;
  end
  Z=d(:,Px0:Px1,:);
  switch strcat('Z',num2str(i))        %ƴ��Z1~Z7
      case 'Z1'
          PIN0=Z;
      case 'Z2'
          PIN1=Z;
      case 'Z3'
          PIN2=Z;
      case 'Z4'
          PIN3=Z;
      case 'Z5'
          PIN4=Z;
      case 'Z6'
          PIN5=Z;
      otherwise 
          PIN6=Z;
  end
  figure(11);
  subplot(1,7,i);
  imshow(Z);
  Px0=Px1;
end
%%   ����ѵ���õ�������
load('-mat','ym');
%%  ���ָ��ĳ����ַ���һ������

for ii=1:7   
    switch ii
      case 1
          PIN0=pretreatment(PIN0);
          Z2=A;
          imwrite(Z2,'A.jpg');
      case 2
         PIN1=pretreatment(PIN1);
         Z2=A;
      case 3
          PIN2=pretreatment(PIN2);
          Z2=A;
      case 4
          PIN3=pretreatment(PIN3);
          Z2=A;
      case 5
          PIN4=pretreatment(PIN4);
          Z2=A;
      case 6
          PIN5=pretreatment(PIN5);
          Z2=A;
      otherwise 
          PIN6=pretreatment(PIN6);
          Z2=A;
    end
    
     figure(12);
  subplot(1,7,ii);
  imshow(Z2);   
end

P0=[PIN0',PIN1',PIN2',PIN3',PIN4',PIN5',PIN6'];







%%  ��ʼʶ���ַ�
for i=2:7     %ѭ������ʶ������������ĸ
  T0= sim(net ,P0(:,i));
  T1 = compet (T0) ;
  d =find(T1 == 1)-1
 if (d==0)
    str='A';
 elseif (d==1)
     str='B';
 elseif (d==2)
     str='C';
 elseif (d==3)
     str='D';
 elseif (d==4)
     str='E';
 elseif (d==5)
     str='F';
 elseif (d==6)
     str='G';
 elseif (d==7)
     str='H';
 elseif (d==8)
     str='I';
 elseif (d==9)
     str='J';
 elseif (d==10)
     str='K';
 elseif (d==11)
     str='L';
 elseif (d==12)
     str='M';
 elseif (d==13)
     str='N';
 elseif (d==14)
     str='O';
 elseif (d==15)
     str='P';
 elseif (d==16)
     str='Q';
 elseif (d==17)
     str='R';
 elseif (d==18)
     str='S';
 elseif (d==19)
     str='T';
 elseif (d==20)
     str='U';
 elseif (d==21)
     str='V';
 elseif (d==22)
     str='W';
 elseif (d==23)
     str='X';
 elseif (d==24)
     str='Y';
 elseif (d==25)
     str='Z';
 else
    str=num2str(d);
 end
 switch i       %������λ��
     case 2
         str1=str;
     case 3
         str2=str;
     case 4
         str3=str;
     case 5
         str4=str;
     case 6
         str5=str;
     otherwise
         str6=str;
  end
end 
s=strcat('��',str1,'.',str2,str3,str4,str5,str6); 
figure(13);
imshow('dw.jpg'),title(s);


