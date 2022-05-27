function numArray = AH_cat2num(catArray)
% This function will convert categorical array into numerical array based
% on the label of data (not the number of category)
% AH 2020/5/28

labels = categories(catArray);
idArray = double(catArray);
celArray = labels(idArray);
numArray = cellfun(@str2num,celArray);
end