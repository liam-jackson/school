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