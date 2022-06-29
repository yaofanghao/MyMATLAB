function Q = calculate_parameter_Q(mu_hat, Up, sigma_hat, xlim)
%   新定义的Q值计算方法

uncertainty_value = Up * sigma_hat;
Q = (mu_hat - uncertainty_value - xlim)/uncertainty_value;

end

