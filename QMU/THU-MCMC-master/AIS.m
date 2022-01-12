function logZZ=AIS(parameter_W,parameter_b,parameter_a,beta)
load test;
%%
b1_A=baserate(testbatchdata);
W=parameter_W;
[v_num,h_num]=size(W);
run_num = 100;
b_A = repmat(b1_A,run_num,1); 
a_B = repmat(parameter_a,run_num,1);
b_B = repmat(parameter_b,run_num,1);
logw = zeros(run_num,1);   %�洢wȡ���������
pv_cur = repmat(sigm(b1_A),run_num,1);  %�˴�v_cur�ǹ�һ����ĸ���
v_cur = pv_cur > rand(run_num,v_num);  %���ɳ�ʼvֵ
logw  =  logw - (v_cur*b1_A' + h_num*log(2));   
%%��ʼ������
Wh = v_cur*W + a_B;
Bb_Base = v_cur*b1_A';
Bv = v_cur*parameter_b';
%%% The CORE of an AIS RUN %%%%%
for bb = beta(2:end-1);
    expWh = exp(bb*Wh);
    logw  =  logw + (1-bb)*Bb_Base + bb*Bv + sum(log(1+expWh),2);   
    poshidprobs = expWh./(1 + expWh);    %���p(h)
    poshidstates = poshidprobs > rand(run_num,h_num);
    %��P(h)����p(v)
    pv_cur = 1./(1 + exp(-(1-bb)*b_A - bb*(poshidstates*W' + b_B)));
    v_cur = pv_cur > rand(run_num,v_num);
    %��������
    Wh      = v_cur*W + a_B;
    Bb_Base = v_cur*b1_A';
    Bv      = v_cur*parameter_b';   
    expWh = exp(bb*Wh);
    logw  =  logw - ((1-bb)*Bb_Base + bb*Bv + sum(log(1+expWh),2));  
end
expWh = exp(Wh);
logw  = logw +  v_cur*parameter_b' + sum(log(1+expWh),2);
r_AIS = logsum(logw(:)) -  log(run_num);
logZZ_base = sum(log(1+exp(b1_A))) + (h_num)*log(2);
logZZ = r_AIS + logZZ_base;
