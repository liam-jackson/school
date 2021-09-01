% clear all, close all, clc;

data_filename = 'student_data.csv';
data_table = rows2vars(readtable(data_filename,'ReadRowNames',true));

hours = data_table.Hours;
outcome = data_table.Pass; 
N = size(data_table, 1); % number of students

%%% Part 1. Plotting the contour map of the cross-entropy cost fxn J
w0max = 1;
w1max = 1;

w0_range = -w0max:0.5:w0max;
w1_range = -w1max:0.5:w1max;
[W0, W1] = meshgrid(w0_range, w1_range);

J = zeros(size(W0));
b = 0;
for w0_ind = 1:length(w0_range)
    for w1_ind = 1:length(w1_range)
        w0_temp = w0_range(w0_ind);
        w1_temp = w1_range(w1_ind);
        J_temp = 0;
        for i = 1:N
            x_temp = hours(i);
            y_temp = outcome(i);
            fx = 1 / (1 + exp(-(w0_temp + w1_temp * x_temp)));
%             disp(fx)
            b = [b, fx];
            J(w0_ind, w1_ind) = (-1/N) * (J_temp + y_temp*log(fx) + (1-y_temp)*log(1-fx));
        end        
    end
end