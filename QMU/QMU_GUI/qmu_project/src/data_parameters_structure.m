classdef data_parameters_structure
    properties
        origin_data;
        
        Distribution;
        Distribution_value;
        
        Stimulus_Xi;
        Ni_Stimulus_Xi;
        stage_mitrix_i;
        explode_or_unexplode_data_flag;
        Stimulus_median_X0;
        parameter_Step_sized;
        
        parameter_n;
        parameter_A;
        parameter_B;
        parameter_M;
        parameter_little_b;
        data_validity;
        
        prob_p;
        Up;
        
        parameter_miu;
        parameter_rou;
        parameter_G;
        parameter_H;
        
        Std_Dev_sigma;
        sigma_mean;
        sigma_variance;
        
        Confidence_level;
        
        design_safety_X0;
        explode_prob_Xp;
        sigma_explode_prob_Xp;
        explode_prob_Xp_QuantileU;
        explode_prob_Xp_QuantileU_lowerL;
        Qmu;
        Qmu_M;
        Qmu_U;
        Qmu_Q;
        explode_prob_Xp_upper;
        explode_prob_Xp_lower;
        explode_prob_Xp_QuantileU_upper;
        explode_prob_Xp_QuantileU_lower;
    end
    methods
        % 构造函数 把原始数据赋给 origin_data 然后 确定分布
        % 然后直接计算不用查表的值
        % 确定 n A B M b data_validity
        function obj = data_parameters_structure(data_to_process,n)
            obj.origin_data = data_to_process;
            if ( n == 1)
                obj.Distribution = 'Gaussian distribution or Logarithmic Gaussian distribution';
                obj.Distribution_value = 1;
            else
                obj.Distribution = 'Logistic distribution or Logarithmic Logistic distribution';
                obj.Distribution_value = 0;
            end
            % 排序
            data_to_process = sortrows(obj.origin_data,1);
            obj.Stimulus_Xi = data_to_process(:,1);
            
            [obj.Ni_Stimulus_Xi,obj.explode_or_unexplode_data_flag] = ...
                get_Ni_Stimulus_Xi(data_to_process);
            [obj.stage_mitrix_i,obj.Stimulus_median_X0,obj.parameter_Step_sized] =...
                get_stage_mitrix_and_X0(obj.Stimulus_Xi);
            
            nABMb = calculate_parameter_n_A_B_M_b(obj.Ni_Stimulus_Xi,obj.stage_mitrix_i);
            obj.parameter_n = nABMb.parameter_n;
            obj.parameter_A = nABMb.parameter_A;
            obj.parameter_B = nABMb.parameter_B;
            obj.parameter_M = nABMb.parameter_M;
            obj.parameter_little_b = nABMb.parameter_little_b;
            obj.data_validity = ...
                check_for_data_validity(obj.Stimulus_Xi,...
                obj.parameter_M,obj.Distribution_value);
            
            obj.parameter_miu = calculate_parameter_miu(...
                obj.Stimulus_median_X0,obj.explode_or_unexplode_data_flag,...
                obj.parameter_A, obj.parameter_n, obj.parameter_Step_sized);
        end
        %% 外界查表 来确定parameter_G 和 parameter_H
        %   设置外部值
        function obj = set_the_rou_G_H(obj,parameter_rou,parameter_G,parameter_H)
            obj.parameter_rou = parameter_rou;
            obj.parameter_G = parameter_G;
            obj.parameter_H = parameter_H;
        end
        function obj = set_the_prob_p(obj,prob_p)
            obj.prob_p = prob_p;
            obj.Up = abs(get_Up(prob_p));
        end
        function obj = set_Confidence_level(obj,Confidence_level)
            obj.Confidence_level = Confidence_level;
        end
        function obj = set_design_safety_X0(obj,design_safety_X0)
            obj.design_safety_X0 = design_safety_X0;
        end
        %%
        % 通过查表值来
        function obj = ...
                set_orther_parameters(obj)
            
            obj.Std_Dev_sigma = calculate_Std_Dev_sigma(...
                obj.parameter_rou,obj.parameter_Step_sized);
            
            [obj.sigma_mean,obj.sigma_variance] = ...
                calculate_sigma_mean_variance(obj.parameter_G,...
                obj.parameter_H,obj.parameter_n,obj.Std_Dev_sigma);
            obj.explode_prob_Xp = ...
                calculate_explode_prob_Xp(obj.parameter_miu,obj.prob_p,obj.Std_Dev_sigma);
            obj.sigma_explode_prob_Xp = ...
                calculate_sigma_explode_prob_Xp(obj.sigma_mean,obj.sigma_variance,obj.prob_p);
            obj.explode_prob_Xp_QuantileU = ...
                calculate_explode_prob_Xp_QuantileU(obj.explode_prob_Xp(2),...
                obj.Confidence_level,obj.sigma_explode_prob_Xp);
            
            obj.explode_prob_Xp_QuantileU_lowerL = ...
                calculate_explode_prob_Xp_QuantileU(obj.explode_prob_Xp(1),...
                obj.Confidence_level,obj.sigma_explode_prob_Xp);
            
            obj.Qmu = get_Qmu_struct...
                (obj.explode_prob_Xp_QuantileU,obj.explode_prob_Xp,obj.design_safety_X0);
            obj.Qmu_Q = obj.Qmu.Qmu_Q;
            obj.Qmu_M = obj.Qmu.Qmu_M;
            obj.Qmu_U = obj.Qmu.Qmu_U;
            obj.explode_prob_Xp_upper = obj.explode_prob_Xp(2);
            obj.explode_prob_Xp_lower = obj.explode_prob_Xp(1);
            obj.explode_prob_Xp_QuantileU_upper = obj.explode_prob_Xp_QuantileU(2);
            obj.explode_prob_Xp_QuantileU_lower = obj.explode_prob_Xp_QuantileU(1);
        end
        function p = get_plot(obj)
            x_miritx = [obj.explode_prob_Xp_upper,obj.explode_prob_Xp_QuantileU_upper,obj.design_safety_X0];
            %x_miritx = [Xp_upper,Xpu_upper,X0];
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
            rymax = 1.01*max(x_miritx);
            rymin = 0.95*min(x_miritx);
            axis([0 11 rymin rymax]);
        end
    end
end