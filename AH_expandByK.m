function arrayOut = AH_expandByK(arrayIn,k,bound)
% this function will expand the arrayIn by neighbering k, and find unique
% elements
% if k=5, front expand 5, back expand 2*k=10
a = 1; %front multiplier
b = 1; %back multipler
if nargin < 2
    k = 15; %=0.0005*30000
end
if nargin < 3
    bound = [];
end
if isa(arrayIn, 'gpuArray') % if gpuArray, arrayAll has to be gpuArray too
    arrayAll = gpuArray(ones(numel(arrayIn)*((a+b)*k+1),1));
else
    arrayAll = ones(numel(arrayIn)*((a+b)*k+1),1);
end
for i = 1:numel(arrayIn)
    element = arrayIn(i);
    arrayAll(((i-1)*((a+b)*k+1)+1):((i-1)*((a+b)*k+1)+(a+b)*k+1)) = [(element-a*k):(element+b*k)];
end
arrayOut = unique(arrayAll);
if ~isempty(bound)
    arrayOut(arrayOut<bound(1) | arrayOut>bound(2)) = []; % delete elements outside bound
end