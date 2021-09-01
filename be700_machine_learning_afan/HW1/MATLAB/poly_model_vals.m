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
