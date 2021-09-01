close all; clear all; clc;

data = readtable('cancer_data.csv');

% Href = shan_entropy(data.(width(data)))
% 
% for col = 1:width(data)-1
%     fprintf('\nthe following entropies are for %s\n', data.Properties.VariableNames{col})
%     Havg = avg_entropy(data.(col), data.class)
%     info_gain = Href - Havg
% end

% [ig, ig_struct] = info_gain(data); 

dig = data(data.shape == "irregular",... 
    find(data.Properties.VariableNames ~= "shape" ));
[ig2, ig_struct2] = info_gain(dig);
disp('breakbreakbreak')

% dig2 = dig(dig.radius == "neutral",...
%     find(dig.Properties.VariableNames ~= "color"))

% [ig3, ig_struct3] = info_gain(dig2);
%%
close all; clear all; clc;

data = readtable('train.csv');

data = data(:,{'Pclass','Sex','Parch', 'Survived'});

% [ig, ig_struct] = info_gain(data);
dig = data(data.Sex == "female",... 
    find(data.Properties.VariableNames ~= "Sex" ));
% [ig2, ig_struct2] = info_gain(dig);
% 
dig2 = dig(dig.Pclass == 3,...
    find(dig.Properties.VariableNames ~= "Pclass"));

[ig3, ig_struct3] = info_gain(dig2);



function H = shan_entropy(feature_or_class)
feat_cat = categorical(feature_or_class);

num_states = length(unique(feat_cat));
[uniq_cnt, uniq_feat] = hist(feat_cat, unique(feat_cat));

total_obs = sum(uniq_cnt);

H = 0;
for i = 1:num_states    
    fprintf('the portion of class %s is %s/%s\n', string(uniq_feat{i}), string(uniq_cnt(i)), string(total_obs));
    p_i = uniq_cnt(i) / total_obs;
    H_temp = - p_i * log2(p_i);
    H = H + H_temp;
end
	
end

function Havg = avg_entropy(feature, classif)

feat_cat = categorical(feature);
class_cat = categorical(classif);

[uniq_cnt, uniq_feat] = hist(feat_cat, unique(feat_cat));

total_obs = length(feat_cat);
num_feats = length(uniq_feat);

Havg = 0;

for f = 1:num_feats
    fprintf('\nThe following data is for state: %s\n', string(uniq_feat{f}));
    feat_idx = feat_cat == uniq_feat(f);
    p_i_avg = uniq_cnt(f) / total_obs;
    shan_temp = shan_entropy(classif(feat_idx))
    Havg = Havg + p_i_avg * shan_temp;
end
end

function [ig, ig_struct] = info_gain(data_table)

class_col = data_table.(width(data_table));
fprintf('The reference entropy value follows:\n');
Href = shan_entropy(class_col)

ig_struct = struct();
for col = 1:width(data_table)-1
    ig_struct_temp = struct();
    var_name = data_table.Properties.VariableNames{col};
    fprintf('\nthe following data is for %s\n', var_name)
    Havg = avg_entropy(data_table.(col), class_col)
    info_gain = Href - Havg
    
    ig_struct_temp.feat_name = var_name;
    ig_struct_temp.info_gain_value = info_gain;
    ig_struct(col).feat_name = ig_struct_temp;
    
end
ig = info_gain;
end
