clear all, close all, clc;

A = readtable('biopsy_data_missing_values.csv', 'NumHeaderLines', 1);

%%% Part A - Data Formatting and Cleaning 

var2_mt = find(strcmp(A.Var2, ''));
for i = 1:length(var2_mt)
    ii = var2_mt(i);
    A.Var2{ii} = 'Irregular';
end

var1_nan = find(~isfinite(A.Var1));
for i = var1_nan
    A.Var1(i) = i;
end 

%%% Part B - Naive Bayes

new_data = table();
new_data.Var1 = [1; 2; 3; 4; 5];
new_data.Var2 = {'Irregular'; 'Irregular'; 'Circle'; 'Circle'; 'Triangle'}; 
new_data.Var3 = {'Large'; 'Small'; 'Large'; 'Large'; 'Large'};
new_data.Var4 = {'Convex'; 'Flat'; 'Concave'; 'Convex'; 'Concave'};
new_data.Var5 = {'Rough'; 'Rough'; 'Smooth'; 'Smooth'; 'Smooth'};
new_data.Var6 = {'Neutral'; 'Red'; 'Neutral'; 'Dark'; 'Neutral'};


mal_inds = find(strcmp(A.Var7, 'Malignant'));
ben_inds = find(strcmp(A.Var7, 'Benign'));

NewBiopProbData = struct();
for biop = 1:height(new_data)
    
    sample = new_data(biop, :);

    BiopProbs = struct();
    temp_mal_probs = [];
    temp_ben_probs = [];
    
    for var_num = 2:length(sample.Properties.VariableNames)
        
        var = string(sample.Properties.VariableNames(var_num));
        biop_var_val = string(table2array(sample(:,var_num)));

        mal_cond_inds = find(strcmp(A.Var7, 'Malignant') & strcmp(A.(var), biop_var_val));
        p_var_giv_mal = length(mal_cond_inds)/length(mal_inds);
        temp_mal_probs = [temp_mal_probs, p_var_giv_mal];
        BiopProbs.(strcat('P_of_', biop_var_val, '_given_Mal')) = p_var_giv_mal;
        
        ben_cond_inds = find(strcmp(A.Var7, 'Benign') & strcmp(A.(var), biop_var_val));
        p_var_giv_ben = length(ben_cond_inds)/length(ben_inds);
        temp_ben_probs = [temp_ben_probs, p_var_giv_ben];
        BiopProbs.(strcat('P_of_', biop_var_val, '_given_Ben')) = p_var_giv_ben;

    end
    p_mal = length(mal_inds)/height(A);
    temp_mal_probs = [temp_mal_probs, p_mal];
    
    p_ben = length(ben_inds)/height(A);    
    temp_ben_probs = [temp_ben_probs, p_ben];
    
    BiopProbs.('P_big_pos_Mal') = prod(temp_mal_probs);
    BiopProbs.('P_big_neg_Ben') = prod(temp_ben_probs);
    BiopProbs.('log_P_big_pos_Mal') = log(prod(temp_mal_probs));
    BiopProbs.('log_P_big_neg_Ben') = log(prod(temp_ben_probs));
    
    if BiopProbs.('P_big_pos_Mal') > BiopProbs.('P_big_neg_Ben')
        BiopProbs.('Predicted_Class') = {'Malignant'};
    elseif BiopProbs.('P_big_pos_Mal') < BiopProbs.('P_big_neg_Ben')
        BiopProbs.('Predicted_Class') = {'Benign'};
    else
        BiopProbs.('Predicted_Class') = {'Inconclusive'};
    end
    
    NewBiopProbData.(['biop' num2str(biop)]) = BiopProbs;
    clear BiopProbs temp_mal_probs temp_ben_probs;
    
end

fields = fieldnames(NewBiopProbData);
for biop_num = 1:length(fields)
    label = ['biop' num2str(biop_num)];
    biop_prob_set = NewBiopProbData.(['biop' num2str(biop_num)]);
    
    disp(label)
    disp(biop_prob_set)
end

