%% Part 1
close all; clear all; clc;

A = imread('default_rgb_reference.tif');
B = double(A);

r = B(:,:,1);
r_flat = reshape(r, [], 1);
g = B(:,:,2);
g_flat = reshape(g, [], 1);
b = B(:,:,3);
b_flat = reshape(b, [], 1);

figure();
plot3(r_flat, g_flat, b_flat, 'k.','MarkerSize', 20);
title('RGB Intensity distribution for pixels in default pic');
xlabel('R');
ylabel('G');
zlabel('B');
axis([0, 255, 0, 255, 0, 255]);
grid on;
hold off;

k = 3; 
X = [r_flat, g_flat, b_flat];
[idx, cent] = kmeans(X, k, 'Replicates', 10);

figure();
cluster_data = struct('cluster',{}, 'rgb_data', {});
col = {'red', 'green', 'blue'};
for clust_num = 1:k
    cluster_data(clust_num).cluster = clust_num;
    rgb_data_temp = X(find(idx == clust_num),:);
    cluster_data(clust_num).rgb_data = [find(idx == clust_num), rgb_data_temp];
    plot3(rgb_data_temp(:,1), rgb_data_temp(:,2), rgb_data_temp(:,3),...
        '.', 'Color', string(col(clust_num)), 'MarkerSize', 20);
    hold on;
end
title(sprintf('(k = %d)-means clustering on default pic', k));
xlabel('R');
ylabel('G');
zlabel('B');
axis([0, 255, 0, 255, 0, 255]);
grid on;
hold off;
clear clust_num;

figure();

for clust_num = 1:k
    im_data = NaN(size(A));
    r_layer = im_data(:,:,1);
    g_layer = im_data(:,:,2);
    b_layer = im_data(:,:,3);
    
    clust_idx = cluster_data(clust_num).rgb_data(:,1);
    r_layer(clust_idx) = r_flat(clust_idx);
    g_layer(clust_idx) = g_flat(clust_idx);
    b_layer(clust_idx) = b_flat(clust_idx);
    
    im_data = cat(3, r_layer, g_layer, b_layer);
    subplot(k, 1, clust_num); 
    image(uint8(im_data));
    title(sprintf('Filtered layer for cluster %d',clust_num)); 
end

%% Part 2

close all; clear all; clc;

%%% Part 1
A = imread('confocal_image01.tif');
B = double(A);

r = B(:,:,1);
r_flat = reshape(r, [], 1);
g = B(:,:,2);
g_flat = reshape(g, [], 1);
b = B(:,:,3);
b_flat = reshape(b, [], 1);

k = 4; 
X = [r_flat, g_flat, b_flat];
[idx, cent] = kmeans(X, k, 'Replicates', 50);

figure();
cluster_data = struct('cluster',{}, 'rgb_data', {});
col = {'red', 'green', 'blue', 'magenta', 'cyan'};
for clust_num = 1:k
    cluster_data(clust_num).cluster = clust_num;
    rgb_data_temp = X(find(idx == clust_num),:);
    cluster_data(clust_num).rgb_data = [find(idx == clust_num), rgb_data_temp];
    plot3(rgb_data_temp(:,1), rgb_data_temp(:,2), rgb_data_temp(:,3),...
        '.', 'Color', string(col(clust_num)), 'MarkerSize', 20);
    hold on;
end
title(sprintf('(k = %d)-means clustering on confocal pic', k));
xlabel('R');
ylabel('G');
zlabel('B');
axis([0, 255, 0, 255, 0, 255]);
grid on;
hold off;
clear clust_num;

figure();

for clust_num = 1:k
    im_data = NaN(size(A));
    r_layer = im_data(:,:,1);
    g_layer = im_data(:,:,2);
    b_layer = im_data(:,:,3);
    
    clust_idx = cluster_data(clust_num).rgb_data(:,1);
    r_layer(clust_idx) = r_flat(clust_idx);
    g_layer(clust_idx) = g_flat(clust_idx);
    b_layer(clust_idx) = b_flat(clust_idx);
    
    im_data = cat(3, r_layer, g_layer, b_layer);
    subplot(round(k/2), 2, clust_num); 
    image(uint8(im_data));
    title(sprintf('Filtered layer for cluster %d',clust_num));
end

disp(sprintf('(k = %d) - clusters needed', k))

%% Part 3


close all; clear all; clc;

%%% Part 1
A = imread('Bacteria_image01.tif');
B = double(A);

r = B(:,:,1);
r_flat = reshape(r, [], 1);
g = B(:,:,2);
g_flat = reshape(g, [], 1);
b = B(:,:,3);
b_flat = reshape(b, [], 1);

k = 3; 
X = [r_flat, g_flat, b_flat];
[idx, cent] = kmeans(X, k, 'Replicates', 50);

figure();
cluster_data = struct('cluster',{}, 'rgb_data', {});
col = {'red', 'green', 'blue', 'magenta', 'cyan'};
for clust_num = 1:k
    cluster_data(clust_num).cluster = clust_num;
    rgb_data_temp = X(find(idx == clust_num),:);
    cluster_data(clust_num).rgb_data = [find(idx == clust_num), rgb_data_temp];
    plot3(rgb_data_temp(:,1), rgb_data_temp(:,2), rgb_data_temp(:,3),...
        '.', 'Color', string(col(clust_num)), 'MarkerSize', 20);
    hold on;
end
title(sprintf('(k = %d)-means clustering on bacteria pic', k));
xlabel('R');
ylabel('G');
zlabel('B');
axis([0, 255, 0, 255, 0, 255]);
grid on;
hold off;
clear clust_num;

figure();

for clust_num = 1:k
    im_data = NaN(size(A));
    r_layer = im_data(:,:,1);
    g_layer = im_data(:,:,2);
    b_layer = im_data(:,:,3);
    
    clust_idx = cluster_data(clust_num).rgb_data(:,1);
    r_layer(clust_idx) = r_flat(clust_idx);
    g_layer(clust_idx) = g_flat(clust_idx);
    b_layer(clust_idx) = b_flat(clust_idx);
    
    im_data = cat(3, r_layer, g_layer, b_layer);
    subplot(round(k/2), 2, clust_num); 
    image(uint8(im_data));
    title(sprintf('Filtered layer for cluster %d',clust_num));
    
end

disp(sprintf('(k = %d) - clusters needed', k))
