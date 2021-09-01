feat = struct();
for i = 1:3
    n = ['cir','irr','tri'];
    a = ['a','b','c','d'];
    feat(i).name = n(i);
    
end

data1 = struct('x',[1,2,3],'y',[4,5,6]); %data structure
b = struct('name',{'Bob'},'data',{data1});