%% ѧϰĿ�꣺�������Ʊ�ϵͳģ��������ϵͳ

load('trainData.mat')    %��������
trainData=[X Y];         %10��������Ϊѵ��
in_format=genfis1(trainData);
[format1,error1,stepsize]=anfis(trainData,in_format,100);    %ѵ��100��
y1=evalfis(X,format1)                  %��ѵ����������
y2=evalfis([42 14 12 786 3936],format1)%��������������
plot(1:10,Y,'b--',1:10,y1,'k.')
legend('Y:ѵ������','y1:��֤����' );   
xlabel('�������'),
ylabel('ϵͳ���')  
%%   ����QQ��1960009019
%%   ���߽���΢�Ź��ںţ�����һƷ��