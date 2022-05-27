function out = getFieldSize(S,n)
% The goal of the function is to get the size of a field of a structure
% without knowing the name of the field
% be careful the index is based on the order the field was created, not by
% alphabet
Cell = struct2cell(S);
dat = Cell{n};
out = size(dat);
