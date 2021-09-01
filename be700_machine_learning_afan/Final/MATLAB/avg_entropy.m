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
    Havg = Havg + p_i_avg * shan_entropy(classif(feat_idx));
end
end