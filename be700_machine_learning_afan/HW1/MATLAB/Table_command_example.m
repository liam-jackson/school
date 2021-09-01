close all
clear all;


% --  A quick table example !!  


% -- Suppose we have want echo a table containing some bass fishing
%    stats:  fishID, fish length, and fish weight 
%
%
%    Also, suppose you forgot to measure the fish weight for 2 fish.
%    We will use "NaN"  (not a number) as a space filler for the
%    table !
%
%    Note:  "NaN" is actually a numerical quantity.... and not a string !
%           Therefore, it can be placed inside a matrix with real
%           numbers in it.
%
%
%    What we want to echo is:
%
%
%          FishID                 length                 weight 
% --------------------------------------------------------------------
%            1                     2.6                     5.1
%            2                     4.7                     8.2
%            3                     3.6                     NaN
%            4                     1.6                     NaN
%
%

% -- Create a string matrix for your field names

my_table_fields = {'FishID'
                   'length'
                   'weight'  };
               

% -- Record the fish parameters

fish_ID = [1 2 3 4]';

fish_length = [2.6  4.7  3.6  1.6]';

fish_weight = [5.1  8.2  NaN  NaN]';


% -- Now, create table !

my_fish_table = table(fish_ID, fish_length, fish_weight, 'VariableNames', my_table_fields)


% -- If you have chosen you table fields names ('FishID', 'length',
%    'weight") in a nice way, you can use them in a struct format to
%    echo individual columns within your table !


% -- Example #1:  I want to echo the fish weights only

my_fish_table.weight


% -- Example #2:  We can calculuate the average of the 2 known fish
%                 weights using the "nanmean" function !
%
%                 "nanmean"  =  Will ignore NaN placeholders and calculate
%                               the averages of the rest of the entries

mean_weight =  nanmean(my_fish_table.weight)














