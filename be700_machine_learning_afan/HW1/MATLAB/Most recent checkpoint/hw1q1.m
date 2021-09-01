% Liam Jackson 
% BE700 ML -- HW1Q1

%%% Ok, honestly this code is overcomplicated. I had to rewrite it 5
%%% times because MATLAB didn't save my changes a couple days in a row. So
%%% in the interest of time, I tried writing functions to accomplish the
%%% analysis for question 1 that I could then recycle for question 2. It's
%%% ugly, I switch data structure types all over the place. But hopefully I
%%% approached the correct answers in the end. 

close all, clear all, clc;
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
warning('off','MATLAB:nearlySingularMatrix')

%%%%%%%%%% Part 1 %%%%%%%%%%

%%% Importing / Sorting Data
[x1, x2, y] = textread('besseldata.txt', ' %f%f%f', 'headerlines', 1);
r = sqrt(x1.^2 + x2.^2);
r_norm = normalize(r);

data_arr = sortrows([x1, x2, r, r_norm, y], 3); 
data_table = array2table(data_arr,...
    'VariableNames', {'x1','x2','r','r_norm','y'});

%%% Bessel Approx
k_bes = 1;
bes_approx = besselj(0, k_bes*data_table.r);

fig1 = figure(1);
dot_sz = 0.2;
line_w = 2.5;
scatter(data_table.r, data_table.y, dot_sz, '.');
hold on;
plot(data_table.r, bes_approx, 'LineWidth', line_w);
hold off;
title({'Timpanic Memb Displacement', 'approximated by Bessel Fxn (J_0)'});
xlabel('r');
ylabel('Intensity');
legend({'Real Data', 'J_0'});

%%% Polynomial Approximations
max_poly_order = 14;
model_data = poly_model_vals(data_table, max_poly_order);
ry_polyvals_table = model_data.ry_polyvals_table;

%%% Calculate Residuals
residuals_table = res_table(ry_polyvals_table)

%%% Plotting LS Poly fits
data_labels = ry_polyvals_table.Properties.VariableNames;

fig2 = figure(2);
sgtitle({'Membrane Displacement Data', 'vs. OLS Polynomial Fits'});
number_of_plots = max_poly_order;

for plot_id = 1:number_of_plots
    subplot(number_of_plots / 2, 2, plot_id);
    scatter(ry_polyvals_table.r, ry_polyvals_table.y, dot_sz, '.');
    hold on;
    plot(ry_polyvals_table.r, ry_polyvals_table.(string(data_labels(plot_id + 2))), 'LineWidth', line_w)
    hold off;
    xlabel('r');
    ylabel('Displacement');
    legend({'y real', string(data_labels(plot_id + 2))});
end

%%%%%%%%%% Part 2 %%%%%%%%%%

cv_rounds = 20;
k_cv = 5;

PE_arr = zeros([cv_rounds, max_poly_order]);
MSE_arr = zeros([max_poly_order, k_cv, cv_rounds]);
all_cv_poly_coeffs = zeros([max_poly_order + 1, max_poly_order, k_cv, cv_rounds]);

for cv_round = 1:cv_rounds 
    binned_data_struct = bin_this_data(data_arr, k_cv);
    binned_data_cell = binned_data_struct.cell;
    bin_indices = 1:k_cv;
    
    for test_bin = 1:k_cv
        train_bins = bin_indices(1:end ~= test_bin);
        train_data_cell = binned_data_cell(train_bins);
        
        train_data_arr = sortrows(cat(1, train_data_cell{:}), 3);
        test_data_arr = sortrows(cell2mat(binned_data_cell(test_bin)), 3);        
        
        train_data_table = array2table(train_data_arr,...
            'VariableNames', {'x1','x2','r','r_norm','y'});
        test_data_table = array2table(test_data_arr,...
            'VariableNames', {'x1','x2','r','r_norm','y'});
        
        model_train_struct = poly_model_vals(train_data_table, max_poly_order);
        cv_ry_polyvals_table = model_train_struct.ry_polyvals_table;
        cv_coeffs_arr = model_train_struct.coeffs_arr;
        all_cv_poly_coeffs(:, :, test_bin, cv_round) = cv_coeffs_arr;
        
        poly_zeros_pad = zeros([size(test_data_arr, 1), max_poly_order]);
        model_poly_vals = [test_data_arr, poly_zeros_pad];
        for poly_ord_ind = 1:max_poly_order
            n_coeffs = poly_ord_ind + 1;
            temp_coeffs = cv_coeffs_arr(1:n_coeffs, poly_ord_ind);
            model_poly_vals(:, poly_ord_ind + 5) = polyval(temp_coeffs, model_poly_vals(:, 4)); 
        end
        
        temp_ry_polyvals_table = array2table([model_poly_vals(:, 3), model_poly_vals(:, 5), model_poly_vals(:,6:end)],...
            'VariableNames', cv_ry_polyvals_table.Properties.VariableNames);
        temp_residuals_table = res_table(temp_ry_polyvals_table);
        MSE_arr(:, test_bin, cv_round) = temp_residuals_table.MSE; 
        
    end
    PE_col = mean(squeeze(MSE_arr(:, :, cv_round)), 2);
    PE_arr(cv_round, :) = PE_col';
end

PE_var_labels = cv_ry_polyvals_table.Properties.VariableNames;
PE_var_labels = PE_var_labels(3:end);
PE_row_nums = 1:cv_rounds;
PE_row_labels = "rnd" + PE_row_nums;

PE_table = array2table(PE_arr,...
    'VariableNames', PE_var_labels,...
    'RowNames', PE_row_labels)

fig3 = figure(3);
plot(PE_table{:,:}.'); 
title('PE values for 20 rounds of (k=5)-CV');
xlabel('Polynomial Model Order');
ylabel('Predictive Error');
legend(PE_table.Properties.RowNames);

%%%%%%%%%% Part 3 %%%%%%%%%%

% A polynomial OLS-fit of order 10 seems to have the best compromise of
% accuracy and economy of variables. A substantial reduction in error
% occurs from order-9 to order-10, with no substantial decrease with
% additional (11, 12, 13, 14) order terms.

char({'A polynomial OLS-fit of order 10 seems to have the best',...
    'compromise of accuracy and economy of variables. A substantial',...
    'reduction in error occurs from order-9 to order-10, with no',...
    'substantial decrease with additional (11, 12, 13, 14) order terms.'})

%%%%%%%%%% Part 4 %%%%%%%%%%

data_table_opt = data_table;

x1 = data_table_opt.x1;
x2 = data_table_opt.x2;
r = data_table_opt.r;
y = data_table_opt.y;

opt_ord = 10;
beta_deg10 = ols_coeffs(r, y, opt_ord);
beta_ud = flipud(beta_deg10);

[X1, X2] = meshgrid(-20:.2:20);
Y_opt = polyval(beta_ud, sqrt(X1.^2 + X2.^2));
J0 = besselj(k_bes, sqrt(X1.^2 + X2.^2));

fig3 = figure(3);
scatter3(x1, x2, y, 15, '.');
hold on;
opt = surf(X1, X2, Y_opt, 'EdgeColor', 'none');
hold off;
title({'Polynomial (p=10) Model', 'vs. Real Displacement Data'});
xlabel('x1');
ylabel('x2');
zlabel('Displacement');
zlim([-.8, 1.5]);

fig4 = figure(4);
scatter3(x1, x2, y, 15, '.');
hold on;
surf(X1, X2, J0, 'EdgeColor', 'none');
title({'Bessel Fxn J_0', 'vs. Real Displacement Data'});
xlabel('x1');
ylabel('x2');
zlabel('Displacement');

%%%%%%%% Functions %%%%%%%%%

function res_squared = residuals(y_real, y_model)
res_squared = sum(abs(y_real - y_model).^2); 
end

function residuals_table = res_table(ry_polyvals_table)

max_poly_order = size(ry_polyvals_table, 2) - 2;
residuals_arr = zeros([max_poly_order, 3]);

for order_ind = 1:max_poly_order
    y_real = ry_polyvals_table.y;
    y_model = table2array(ry_polyvals_table(:, order_ind + 2));
    residuals_arr(order_ind, 1) = order_ind;
    residuals_arr(order_ind, 2) = residuals(y_real, y_model);
    residuals_arr(order_ind, 3) = residuals_arr(order_ind, 2) ./ size(ry_polyvals_table, 1);
end

residuals_table = array2table(residuals_arr,...
    'VariableNames',{'Polynomial_Order', 'Residual Sum', 'MSE'});
end

function beta = ols_coeffs(r, y, poly_order)

X = zeros(length(r), poly_order + 1);
X(:, 1) = 1;

for ord_ind = 1:poly_order
    X(:, ord_ind + 1) = r.^(ord_ind);
end

beta = (X' * X) \ (X' * y);

end

function model_data_struct = poly_model_vals(data_table, max_poly_order)
model_data_struct = struct();

r = data_table.r;
r_norm = data_table.r_norm;
y = data_table.y;

n_coeffs = max_poly_order + 1;

coeffs_arr = zeros([n_coeffs, max_poly_order]); %15x14
poly_vals = zeros([length(r), max_poly_order]); %5000x14

for poly_order_ind = 1:max_poly_order
    beta = ols_coeffs(r_norm, y, poly_order_ind);    
    poly_coeffs = flipud(beta);
    for r_ind = 1:length(poly_coeffs)
        coeffs_arr(r_ind, poly_order_ind) = poly_coeffs(r_ind); %Array is in DESCENDING ORDER of poly coeffs 
    end
    poly_vals(:, poly_order_ind) = polyval(poly_coeffs, r_norm);    
end

data_poly_vals_arr = [r, y, poly_vals];

data_var_names = {'r', 'y'};
poly_var_nums = 1:length(poly_vals(1,:));
poly_var_names = "p=" + poly_var_nums;
var_names = {[data_var_names, poly_var_names]};

data_poly_vals_table = array2table(data_poly_vals_arr,...
    'VariableNames',var_names{1});

model_data_struct.coeffs_arr = coeffs_arr;
model_data_struct.ry_polyvals_arr = data_poly_vals_arr;
model_data_struct.ry_polyvals_table = data_poly_vals_table;

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
