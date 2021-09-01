clear all; close all; clc;

%%% Part 1
opts = detectImportOptions('housing.csv', 'NumHeaderLines', 1);
preview('housing.csv', opts)

A = readtable('housing.csv', 'NumHeaderLines', 1);
B = A(:,1:end-1);
labels = {'<1H OCEAN', 'INLAND', 'NEAR OCEAN', 'NEAR BAY', 'ISLAND'};

for label_ind = 1:numel(labels)
    for inst = 1:height(A)
        if string(A.Var10(inst)) == labels(label_ind)
            B.Var10(inst) = label_ind;
        end
    end
end

figure();
pltmtx = plotmatrix(table2array(B));

figure();
hold on;
separated_data.label = struct('title',{},'label_data', {});

for label_ind = 1:numel(labels)
    label = string(labels(label_ind));
    separated_data.label(label_ind).title = label;

    all_data_for_label = B(B.Var10 == label_ind,:);
    separated_data.label(label_ind).label_data = all_data_for_label;
    
    scatter(table2array(all_data_for_label(:,1)),...
        table2array(all_data_for_label(:,2)), '.',...
        'DisplayName', label);    
end
title('Locations of 5 Housing Classes across California');
xlabel('Longitude x_1');
ylabel('Latitude x_2');
legend();

%%% Part 2
N = readtable('newhouses.csv', 'NumHeaderLines', 1);

x_new = [N.Var2, N.Var3];
x_old = [B.Var1, B.Var2, B.Var10];

% Overlay New Houses:
scatter(x_new(:,1), x_new(:,2), 1000, 'k.', 'DisplayName', 'New Data');

nn_data.house = struct('title',{},'nns',{}, 'pred', {});

% Find k=20 NNs
for N_idx = 1:height(N)
    nn_idx = knnsearch(x_old(:,1:2), x_new(N_idx,:), 'K', 20); %, 'SortIndices', false);
    nns = x_old(nn_idx,:);
    
    house_name = strcat('House',string(N_idx));
    nn_data.house(N_idx).title = house_name;
    nns_table = array2table(nns, 'VariableNames',{'Longitude','Latitude','Label'});
    nn_data.house(N_idx).nns = nns_table;
    nn_labels = unique(nns(:,end)); 
    label_freq = [nn_labels, histc(nns(:,end),nn_labels)];
    label_pred = label_freq(label_freq(:,2)==max(label_freq(:,2)),1);
    nn_data.house(N_idx).pred = label_pred;
    
    disp(strjoin(['20 NNS of', house_name, 'and their data:'], ' '))
    disp(nns_table)
    
    scatter(nns(:,1),nns(:,2), 500, 'b^', 'DisplayName', 'NNs')
    legend('AutoUpdate','off')
end
axis square;
hold off;

for N_idx = 1:height(N)
    house_name = nn_data.house(N_idx).title;
    label_pred = nn_data.house(N_idx).pred;
    disp(strjoin(['Predicted label for', house_name, 'is', cellstr(labels(label_pred))], ' '))
end

for house_num = 1:height(N)
    wm = webmap('Open Street Map');
    hold on;
    nns_webmap_data = geopoint(nn_data.house(house_num).nns.Latitude, nn_data.house(house_num).nns.Longitude);
    nns_webmarker = wmmarker(nns_webmap_data, 'Color', 'blue');
    new_houses_data = geopoint(x_new(house_num, 2), x_new(house_num, 1));
    new_houses_webmarker = wmmarker(new_houses_data, 'Color', 'red');
end
hold off;










