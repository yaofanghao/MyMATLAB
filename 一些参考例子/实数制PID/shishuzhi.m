%%  GA(Generic Algorithm) Program to optimize PID Parameters
clc
clear all
close all
%%  ����ȫ�ֱ���rin yout timef���⼸����������������Ӻ����ж����Բ���
global rin yout timef

%% ��ֵ
Size=30;
CodeL=3;

%%  ��ʼ��

%%  MinX()������һ��1*3��һά������ֵ��[0,0,0]
MinX(1)=zeros(1);
MinX(2)=zeros(1);
MinX(3)=zeros(1);
%%  MaxX()������һ��1*3��һά������ֵ��[20,1,1]
MaxX(1)=20*ones(1);
MaxX(2)=1.0*ones(1);
MaxX(3)=1.0*ones(1);
%%  KPID������һ��30��3�еľ���
KPID(:,1)=MinX(1)+(MaxX(1)-MinX(1))*rand(Size,1);%��KPID����ĵ�1�и�ֵ
KPID(:,2)=MinX(2)+(MaxX(2)-MinX(2))*rand(Size,1);%��KPID����ĵ�2�и�ֵ
KPID(:,3)=MinX(3)+(MaxX(3)-MinX(3))*rand(Size,1);%��KPID����ĵ�3�и�ֵ

G=100;
BsJ=0;

for kg=1:1:G        %forѭ����kg��1��100��ÿѭ��һ�ξͼ�1
    time(kg)=kg;%ѭ����time�����ÿ��ֵ��ֵ
    
    %ÿ��ѭ��һ��kg�����ﶼҪ���һ��������i��1��Size��ѭ��
    for i=1:1:Size
        KPIDi=KPID(i,:);

        [KPIDi,BsJ]=chap5_3f(KPIDi,BsJ);

        BsJi(i)=BsJ;
    end
    %sort�������Ŵ�С������1�ĵ�һ������ά����������BsJi��Ԫ�ء�
    [OderJi,IndexJi]=sort(BsJi);
    BestJ(kg)=OderJi(1);
    BJ=BestJ(kg);
    Ji=BsJi+1e-10;    
    fi=1./Ji;
    %  Cm=max(Ji);
    %  fi=Cm-Ji;                     
    %sort�������Ŵ�С������1�ĵ�һ������ά����������fi��Ԫ�ء�
    [Oderfi,Indexfi]=sort(fi);       
    Bestfi=Oderfi(Size);          
    BestS=KPID(Indexfi(Size),:);  
    %��һ��ֻ����ʾfi�����е����ֵ���Գ���û��Ӱ��
    max(fi)
   
%    kg   
%    BJ
%    BestS
    %�Ծ���fi���
    fi_sum=sum(fi);
    fi_Size=(Oderfi/fi_sum)*Size;
    %floor������A��Ԫ����������ΪС�ڻ����A�����������
    fi_S=floor(fi_Size);                    
    r=Size-sum(fi_S);
   
    Rest=fi_Size-fi_S;
    %sort�������Ŵ�С������1�ĵ�һ������ά����������BsJi��Ԫ�ء�
    [RestValue,Index]=sort(Rest);
   
    for i=Size:-1:Size-r+1
        fi_S(Index(i))=fi_S(Index(i))+1;     
    end

    k=1;
    for i=Size:-1:1         
        for j=1:1:fi_S(i)  
            TempE(k,:)=KPID(Indexfi(i),:);       
            k=k+1;                            
        end
    end
   
    Pc=0.90;
    for i=1:2:(Size-1)
          temp=rand;
          %���Pc����temp,�ͽ������µĸ�ֵ
      if Pc>temp                      
          alfa=rand;
          TempE(i,:)=alfa*KPID(i+1,:)+(1-alfa)*KPID(i,:);  
          TempE(i+1,:)=alfa*KPID(i,:)+(1-alfa)*KPID(i+1,:);
      end
    end
    TempE(Size,:)=BestS;
    KPID=TempE;
    
    Pm=0.10-[1:1:Size]*(0.01)/Size;      
    Pm_rand=rand(Size,CodeL);
    Mean=(MaxX + MinX)/2; 
    Dif=(MaxX-MinX);

   for i=1:1:Size
      for j=1:1:CodeL
         if Pm(i)>Pm_rand(i,j)        
            TempE(i,j)=Mean(j)+Dif(j)*(rand-0.5);
         end
      end
   end
    %Guarantee TempE(Size,:) belong to the best individual
    %��֤TempE�����е���ֵ�����Ż���õ�
   TempE(Size,:)=BestS;      
   KPID=TempE;
end
%����Ҳֻ�ǽ���ֵ��ӡ�������д��ڣ��Գ���û��ʲôӰ��
Bestfi
BestS
Best_J=BestJ(G);
%�½�һ����ͼ����
figure(1);
%��timeΪ�����꣬BestJΪ�������ͼ
plot(time,BestJ);
%��ͼ�ĺ��������������������
xlabel('Times');ylabel('Best J');
figure(2);%�½�һ����ͼ����figure(2)
%��timefΪ�����꣬rinΪ�������ͼ����ɫΪred��������timefΪ�����꣬youtΪ�������ͼ��������ɫΪblue
plot(timef,rin,'r',timef,yout,'b');
%��ͼ�ĺ��������������������
xlabel('Time(s)');ylabel('rin,yout');
