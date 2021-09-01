close all, clear all, clc;
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

data_table = readtable('winequality_red.csv');

t = data_table.('citricAcid');
y = data_table.('fixedAcidity');

u = (t - mean(t))./std(t); %standardizing t

var_names = {'t','u','y'};
data_table = sortrows(array2table([t, u, y],...
    'VariableNames', var_names), 't');
data_arr = table2array(data_table);

% plotmatrix(data_arr);

p = 9;
k_cv = 5;
alpha_range = 0:50:2500;
num_PE_curves = 10;

PE_arr = [alpha_range', zeros([length(alpha_range), num_PE_curves])]
%%

for PE_curve = 1:1:num_PE_curves
    
    binned_data_struct = bin_this_data(data_arr, k_cv);
    
    %%% celldisp(binned_data_struct.table.bin1) %%% does not display
    %%% ordered array
    
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
%             X = [train_data_table.u];
            D = x2fx(X, 'interaction');
            D(:,1) = []; % No constant term
            
            disp(alpha)
            B = ridge(train_data_table.y, D, alpha)
            beta_flip = flipud(B);
            
            u_test = test_data_table.u;
            y_real = test_data_table.y;
%             y_model = polyval(beta_flip, u_test);
            y_model = polyval(B, u_test);            
            

            
            res_sq_sum = residuals(y_real, y_model);
            MSE_5fold(test_bin, 1) = res_sq_sum / length(y_real); 
            
%             scatter(test_table.u, test_table.yr); hold on; scatter(test_table.u, test_table.ym);

        end
        disp(alpha_ind)
        PE_arr(alpha_ind, PE_curve + 1) = mean(MSE_5fold) 
    end
end
%%

%%% Plot all PE Curves
close all;

for PE_curve_id = 1:num_PE_curves
    plot(PE_arr(:, 1), PE_arr(:, PE_curve_id + 1));
    hold on;
end
    
%%
            test_table = table();
            test_table.u = u_test;
            test_table.yr = y_real;
            test_table.ym = y_model;















