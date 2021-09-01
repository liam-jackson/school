close all, clear all, clc;
warning('off','MATLAB:table:ModifiedAndSavedVarnames')

%%%%%%%%%% Part 1 %%%%%%%%%%

data_table_raw = readtable('winequality_red.csv');

t = data_table_raw.('citricAcid');
u = (t - mean(t))./std(t); %standardizing t
y = data_table_raw.('fixedAcidity');

var_names = {'t','u','y'};
data_table = sortrows(array2table([t, u, y],...
    'VariableNames', var_names), 't');
data_arr = table2array(data_table);

A = table2array(data_table_raw);
fig1 = figure(1);
plotmatrix(A);
title('Scatter Plot Matrix of Wine Data')

%%%%%%%%%% Part 2 %%%%%%%%%%

p = 9;
k_cv = 5;
alpha_range = 0:50:3000;
num_PE_curves = 100;

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
            
            u_train = train_data_table.u;
            y_train = train_data_table.y;
            u_test = test_data_table.u;
            y_test = test_data_table.y;

            X = ols_coeffs_data(u_train, y_train, p).X;
            X(:,1) = [];
            
            B = ridge(y_train, X, alpha, 0);
            beta_flip = flipud(B);
            
            
            y_model = polyval(beta_flip, u_test);
            
            res_sq_sum = residuals(y_test, y_model);
            MSE_5fold(test_bin, 1) = res_sq_sum / length(y_test); 
            
        end
        PE_arr(alpha_ind, PE_curve + 1) = mean(MSE_5fold); 
    end
end

%%% Plot all PE Curves
fig2 = figure(2);
for PE_curve_id = 1:num_PE_curves
    plot(PE_arr(:, 1), PE_arr(:, PE_curve_id + 1));
    hold on;
end
title({'100 Curves of Predictive Error', '(of a 9th order polynomial model)', 'vs. \alpha in L2-Regularization'});    
ylim([0, 10])
xlabel('\alpha value')
ylabel('PE')
hold off;

opt_alpha = 400;
char('It looks the optimal alpha is maybe ~400 ish')

%%%%%%%%%% Part 3 %%%%%%%%%%

u_full = data_table.u;
y_full = data_table.y;

X_full_std = ols_coeffs_data(u_full, y_full, 9).X;
X_full_std(:,1) = [];

w_opt = ridge(y_full, X_full_std, opt_alpha, 0)
w_opt_flip = flipud(w_opt);

u_linear = linspace(min(u_full), max(u_full), 1500);
y_model_opt = polyval(w_opt_flip, u_linear);

w_ols = ols_coeffs(u_full, y_full, 9)
w_ols_flip = flipud(w_ols);

y_model_ols = polyval(w_ols_flip, u_linear);

%%%%%%%%%% Part 4 %%%%%%%%%%

fig3 = figure(3);
scatter(u_full, y_full, 100, '.')
hold on;
plot(u_linear, y_model_ols, 'r--')
plot(u_linear, y_model_opt, 'g-.')
title({'Raw Data', 'vs OLS fit', 'vs L2 Regularized Fit (\alpha = 400)'});
xlabel('u (standardized t) [Citric Acid Content]')
ylabel('Fixed acidity')
legend({'Raw Data', 'OLS (p = 9)', 'L2 Reg'}, 'location', 'southwest')

coeff_table_all = array2table([w_ols, w_opt],...
    'VariableNames', {'w_ols', 'w_opt_L2'},...
    'RowNames', {'p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10'})

%%%%%%%% Functions %%%%%%%%%

function res_squared = residuals(y_real, y_model)
res_squared = sum(abs(y_real - y_model).^2); 
end

function beta = ols_coeffs(r, y, poly_order)

X = zeros(length(r), poly_order + 1);
X(:, 1) = 1;

for ord_ind = 1:poly_order
    X(:, ord_ind + 1) = r.^(ord_ind);
end

beta = (X' * X) \ (X' * y);

end

function ols_data = ols_coeffs_data(x, y, poly_order)
ols_data = struct();

X = zeros(length(x), poly_order + 1);
X(:, 1) = 1;

for ord_ind = 1:poly_order
    X(:, ord_ind + 1) = x.^(ord_ind);
end

beta = (X' * X) \ (X' * y);

ols_data.beta = beta;
ols_data.X = X;
ols_data.res_squares_sum = norm(y - X*beta).^2;

end

function binned_data_struct = bin_this_data(data_arr_to_bin, k_bins)
binned_data_struct = struct();

num_data_pts = size(data_arr_to_bin, 1);
rows_per_bin = floor(num_data_pts / k_bins);
extra_rows_needed = mod(num_data_pts, k_bins);

perm_ind = randperm(num_data_pts);
perm_data = data_arr_to_bin(perm_ind, :);

rowDist = rows_per_bin * ones(1, k_bins);

for extra_row = 1:extra_rows_needed
    rowDist(extra_row) = rowDist(extra_row) + 1;
end

binned_data_cell = mat2cell(perm_data, rowDist)';
bin_nums = 1:k_bins;
bin_names = "bin" + bin_nums;

binned_data_table = cell2table(binned_data_cell,... 
    'VariableNames', bin_names);

binned_data_struct.cell = binned_data_cell;
binned_data_struct.table = binned_data_table;

end









