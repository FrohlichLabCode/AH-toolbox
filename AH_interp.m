function arrayOut = AH_interp(arrayIn)
% This function replace NaN values in a 1d array using its neighboring
% values
% AH 3/9/2021

% First check if there is any missing values
nanIn = isnan(arrayIn);
infIn = isinf(arrayIn);
misIn = nanIn | infIn;
arrayOut = arrayIn;
if sum(misIn) ~= 0 % no missing value
    t = 1:numel(arrayIn);
    arrayOut(misIn) = interp1(t(~misIn), arrayOut(~misIn), t(misIn));
end
% Handle margins
if misIn(1) == 1 % first element is nan
    IDs = find(misIn == 0);
    arrayOut(1:IDs(1)) = arrayOut(IDs(1))*ones(1,IDs(1));
end
if misIn(end) == 1
    IDs = find(misIn == 0);
    arrayOut(end:IDs(end)) = arrayOut(IDs(end))*ones(1,IDs(end));
end