clc;
clear;
data_example = xlsread('./datapaper.xlsx');
dataset = data_parameters_structure(data_example,1);
parameter_rou = 0.52;
parameter_G = 1.134;
parameter_H = 1.146;
dataset = dataset.set_the_rou_G_H(parameter_rou,parameter_G,parameter_H);

dataset = dataset.set_the_prob_p(1e-6);

dataset = dataset.set_Confidence_level(0.9999);

dataset = dataset.set_design_safety_X0(26);

dataset = dataset.set_orther_parameters();

resultTable = generate_result_table(dataset);

p = dataset.get_plot();
