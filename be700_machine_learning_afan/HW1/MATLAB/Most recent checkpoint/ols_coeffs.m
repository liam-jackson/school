function beta = ols_coeffs(r, y, poly_order)

X = zeros(length(r), poly_order + 1);
X(:, 1) = 1;

for ord_ind = 1:poly_order
    X(:, ord_ind + 1) = r.^(ord_ind);
end

beta = (X' * X) \ (X' * y);

end
