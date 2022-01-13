function [meann,stdn] = shengjiangfa(xint, d, mu, sigma)
%SHENGJIANGFA �˴���ʾ�йش˺�����ժҪ
% ģ��������ʵ�飬����xint,d,mu,sigma�����meann,stdn
%   �˴���ʾ��ϸ˵��

n = 50; %������
% xciΪ������̬�����
% ���ʣ�mu��sigma��ֵ���ȷ�������Ƿ�Ӧ����Դ����ʵʵ�����ֵ��
xc = normrnd(mu, sigma, 1, n);
x = zeros(1,n); %���xi
x(1) = xint; %��ʼ�̼���
y = zeros(1,n); %�����Ӧֵ0��1

for i = 1:n
    if x(i) >= xc(i)
        y(i) = 1;
        x(i+1) = x(i)-d;
    else
        y(i) = 0;
        x(i+1) = x(i)+d;
    end
end
x(:,51) = []; %ɾ�����ɵ�x(51)

% �ж�������Ч�ԣ��̼�������ӦΪ4-7��
s = (max(x)-min(x))/d + 1;
if s>=4 || s<=7  
%     fprintf('�̼�������Ϊ%d,������Ч\n', s);
else
    fprintf('�̼�������Ϊ%d,������Ч\n', s); %��ӡ��Ϣ��������Ƿ�����
end

meann = mean(x);
stdn = std(x);

end

