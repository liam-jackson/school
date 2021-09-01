close all, clear all, clc;
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

data_table = readtable('winequality_red.csv');

t = data_table.('citricAcid');
y = data_table.('fixedAcidity');

t_mean = mean(t);
t_sd = std(t);
u = (t - t_mean)./t_sd;

var_names = {'t','u','y'};
data_table = sortrows(array2table([t, u, y],...
    'VariableNames', var_names), 't');
data_arr = table2array(data_table);

% plotmatrix(data_arr);

p = 9;
k_cv = 5;
alpha_range = 0:50:250;
num_PE_curves = 10;

PE_arr = [alpha_range', zeros([length(alpha_range), num_PE_curves])];

for PE_curve = 1:1:num_PE_curves
    
    binned_data_struct = bin_this_data(data_arr, k_cv);
    binned_data_cell = binned_data_struct.cell;
    bin_indices = 1:k_cv;
    
    for alpha = alpha_range   
        alpha_ind = find(alpha_range == alpha);
        MSE_5fold = zeros([k_cv, 1]);
        for test_bin = 1:k_cv
            train_bins = bin_indices(1:end ~= test_bin);
            train_data_cell = binned_data_cell(train_bins);

            train_data_arr = sortrows(cat(1, train_data_cell{:}));
            test_data_arr = sortrows(cell2mat(binned_data_cell(test_bin)));        

            train_data_table = array2table(train_data_arr,...
                'VariableNames', var_names);
            test_data_table = array2table(test_data_arr,...
                'VariableNames', var_names);

            X = ols_coeffs_data(train_data_table.u, train_data_table.y, p).X;
            
            B = ridge(train_data_table.y, X, alpha);
            beta_flip = flipud(B);
            
            u_test = test_data_table.u;
            y_real = test_data_table.y;
            y_model = polyval(beta_flip, u_test);
%             y_model = polyval(B, u_test);            
            
            test_table = table();
            test_table.u = u_test;
            test_table.yr = y_real;
            test_table.ym = y_model;
            
            res_sq_sum = residuals(y_real, y_model)
            MSE_5fold(test_bin) = res_sq_sum / length(u_test)
            scatter(test_table.u, test_table.yr); hold on; scatter(test_table.u, test_table.ym);

        end
        PE_arr(alpha_ind, PE_curve + 1) = mean(MSE_5fold)
    end
end
%%

%%% Plot all PE Curves

for PE_curve_id = 1:size(PE_arr, 2)
    plot(PE_arr(:, 1), PE_arr(:, PE_curve_id+1));
    hold on;
end
    















