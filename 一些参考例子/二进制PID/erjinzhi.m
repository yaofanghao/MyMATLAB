%%  ����ȫ�ֱ���rin yout timef���⼸����������������Ӻ����ж����Բ���
global rin yout timef 
%%   ��������
G=100;     %%  �������ѭ�������й�ϵ��24��
Size=30; %  ������ĵڶ�ʮ�г�������й�ϵ������Ƕ�׵�ѭ��������26��
CodeL=10; %  ������ĵڶ�ʮ�г�������й�ϵ

%%   ��ʼ��

%%   MinX()������һ��1*3��һά������ֵ��[0,0,0]
MinX(1)=zeros(1);
MinX(2)=zeros(1);
MinX(3)=zeros(1);
%%   MaxX()������һ��1*3��һά������ֵ��[20,1,1]
MaxX(1)=20*ones(1);
MaxX(2)=1.0*ones(1);
MaxX(3)=1.0*ones(1);

%%   round��������������������
E=round(rand(Size,3*CodeL));
%%   ��ʼ��BsJ����ֵ��һ���ʼ������0
BsJ=0;
%%  ѭ��
for kg=1:1:G      %ѭ������
    time(kg)=kg;
    for s=1:1:Size   %Ƕ�׵�ѭ������
        m=E(s,:);
        y1=0;y2=0;y3=0;
        m1=m(1:1:CodeL);    %����ǰ������ֵ
        for i=1:1:CodeL     %����ǰ������ֵ
            y1=y1+m1(i)*2^(i-1);
        end
        KPID(s,1)=(MaxX(1)-MinX(1))*y1/1023+MinX(1);
        m2=m(CodeL+1:1:2*CodeL);
        for i=1:1:CodeL
            y2=y2+m2(i)*2^(i-1);
        end
        KPID(s,2)=(MaxX(2)-MinX(2))*y2/1023+MinX(2);
        m3=m(2*CodeL+1:1:3*CodeL);
        for i=1:1:CodeL
            y3=y3+m3(i)*2^(i-1);
        end
            KPID(s,3)=(MaxX(3)-MinX(3))*y3/1023+MinX(3);
            KPIDi=KPID(s,:);
            [KPIDi,BsJ]=chap5_3f(KPIDi,BsJ);
            BsJi(s)=BsJ;
    end

    [OderJi,IndexJi]=sort(BsJi);
    BestJ(kg)=OderJi(1);
    BJ=BestJ(kg);
    Ji=BsJi+1e-10;
    fi=1./Ji;
    [Oderfi,Indexfi]=sort(fi);
    Bestfi=Oderfi(Size);
    BestS=E(Indexfi(Size),:);

    fi_sum=sum(fi);
    fi_Size=(Oderfi/fi_sum)*Size;
    fi_S=floor(fi_Size);
    kk=1;
    for i=1:1:Size
        for j=1:1:fi_S(i)
            TempE(kk,:)=E(Indexfi(i),:);
            kk=kk+1;
        end
    end
    pc=0.60;
    n=ceil(20*rand);
    for i=1:2:(Size-1)
        temp=rand;
        if pc>temp
            for j=n:1:20
            TempE(i,j)=E(i+1,j);
            TempE(i+1,j)=E(i,j);
            end
        end
    end
    TempE(Size,:)=BestS;
    E=TempE;
    pm=0.001-[1:1:Size]*(0.001)/Size; %Bigger fi, smaller pm
    for i=1:1:Size
        for j=1:1:3*CodeL
            temp=rand;
            if pm>temp
                if TempE(i,j)==0
                    TempE(i,j)=1;
                else
                    TempE(i,j)=0;
                end
            end
        end
    end
    TempE(Size,:)=BestS;
 
    E=TempE;
end
%%  ����Ҳֻ�ǽ���ֵ��ӡ�������д��ڣ��Գ���û��ʲôӰ��
Bestfi
BestS
KPIDi
Best_J=BestJ(G)
%%  
figure(1);  %�½�һ����ͼ����
plot(time,BestJ,'r');%��timeΪ�����꣬BestJΪ�������ͼ��������ɫΪred��ɫ
xlabel('Times');ylabel('Best J');%��ͼ�ĺ��������������������
grid on;%��ͼ�μ�������
figure(2);%�½�һ����ͼ����figure(2)
plot(timef,rin,'r',timef,yout,'b');%��timefΪ�����꣬rinΪ�������ͼ����ɫΪ��ɫ��������timefΪ�����꣬youtΪ�������ͼ��������ɫΪblue
xlabel('Time(s)');ylabel('rin,yout');%��ͼ�ĺ��������������������
grid on

