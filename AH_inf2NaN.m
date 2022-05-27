function newTable = AH_inf2NaN(Table)
% converting all Inf into NaN, need to first convert table to
% array then convert back

tmp = Table;
tmp(isinf(tmp)) = NaN;
newTable = tmp;
end