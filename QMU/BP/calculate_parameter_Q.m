function Q = calculate_parameter_Q(mu_hat, Up, sigma_hat, xlim)
%   �¶����Qֵ���㷽��

uncertainty_value = Up * sigma_hat;
Q = (mu_hat - uncertainty_value - xlim)/uncertainty_value;

end

