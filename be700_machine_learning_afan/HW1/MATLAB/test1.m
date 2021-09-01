clear all, clc;


% A = repmat(1:11, 2, 1).'
% n = 4;
% % B = reshape([A(:); zeros(mod(-numel(A),n),1)],n,[])'
% 
% k = 5;
% rowDist = [2 2 2 2 3];
% B = mat2cell(A, rowDist)
% B{:}

data = repmat(1:13, 2, 1).'
k = 5;

num_data_pts = size(data, 1);
rows_per_bin = floor(num_data_pts / k);
extra_rows_needed = mod(num_data_pts, k);

rowDist = rows_per_bin * ones(1, k);

for extra_row = 1:extra_rows_needed
    rowDist(extra_row) = rowDist(extra_row) + 1;
end

perm_data = mat2cell(data, rowDist)
celldisp(perm_data)

