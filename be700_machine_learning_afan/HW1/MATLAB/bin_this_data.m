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





















