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
