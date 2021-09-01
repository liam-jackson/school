function res_squared = residuals(y_real, y_model)
res_squared = sum(abs(y_real - y_model).^2); 
end




