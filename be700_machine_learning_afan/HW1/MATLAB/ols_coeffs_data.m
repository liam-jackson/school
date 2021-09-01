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
