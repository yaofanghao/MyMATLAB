function [parameter_miu, Std_Dev_sigma, explode_prob_Xp, prob, Q] = ...
   shengjiangfa(xint, d, mu, sigma, prob_p, xlim)
%  ģ��������ʵ�飬 ����xint,d,mu,sigma,xlim
%  �������miu-hat��sigma-hat��
%  ����ˮƽΪprob_pʱ����ը������Ӧ��Xp�����䣬��Ӧ����explode_prob_Xp��Q
%  ������1.16

n = 50; %������
% xciΪ������̬�����
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

table = tabulate(x); %��ȡxi����ֵ��Ƶ��
xi = table(:,1); % �̼���xi
sum_v_m = table(:,2); % xi��Ƶ��

%������Ӧ�Ͳ���Ӧ����
v = zeros(length(xi),1);
m = zeros(length(xi),1); 

%ͳ�ƴ̼���Ϊxi(j)ʱ����Ӧ�Ͳ���Ӧ����
for j = 1:length(xi)
    for i = 1:n
        if  x(i) == xi(j)  
            if y(i)==1
                v(j) = v(j) + 1;
            else 
                m(j) = m(j) + 1;
            end
        end
    end
end

% ��ģ�����������ɱ����ʽ
data_example = [xi, v, m, sum_v_m];

Stimulus_Xi = data_example(:,1); %��ȡ�̼���xi

% �Ƚ�v��m��С���ж�nȡvi����mi
[Ni_Stimulus_Xi,explode_or_unexplode_data_flag] = get_Ni_Stimulus_Xi(data_example);

% ��ȡ��λ��X0 �Լ�̨��������,����d
[stage_mitrix_i,Stimulus_median_X0,parameter_Step_sized] =...
    get_stage_mitrix_and_X0(Stimulus_Xi);

% ����A B M b ��ģ��
nABMb = calculate_parameter_n_A_B_M_b(Ni_Stimulus_Xi,stage_mitrix_i);
parameter_n = nABMb.parameter_n;
parameter_A = nABMb.parameter_A;
parameter_B = nABMb.parameter_B;
parameter_M = nABMb.parameter_M;
parameter_little_b = nABMb.parameter_little_b;

% ������Ч���ж� 1��Ч 0 ��Ч
data_validity = check_for_data_validity(parameter_M);
if data_validity == 0
    fprintf('ע�⣡��������Ч����\n');
%     error('MС��0.25�����ɵ�������Ч');
end

% ����mu-hat
parameter_miu = calculate_parameter_miu(...
    Stimulus_median_X0,explode_or_unexplode_data_flag,...
    parameter_A, parameter_n, parameter_Step_sized);

% ����� 
parameter_rou = calculate_parameter_rou(parameter_M,parameter_little_b);
% ����sigma-hat
Std_Dev_sigma = calculate_Std_Dev_sigma(parameter_rou,parameter_Step_sized);

% ��������ˮƽΪprob_pʱ����ը������Ӧ��Xp������ 
Up = abs(get_Up(prob_p));  %��̬��P��λ��
[explode_prob_Xp] = calculate_explode_prob_Xp(parameter_miu,prob_p,Std_Dev_sigma);

%��Ӧ���� ��((mu-hat - mu) / sigma)
prob = cdf('Normal',parameter_miu,mu,sigma);
% prob = (length(find(y==1)))/n;

%�¶����Qֵ���㷽��
Q = calculate_parameter_Q(parameter_miu, Up, Std_Dev_sigma, xlim);

end

