function arrayOut = AH_expandByK2(arrayIn,k,bound)
% this function will expand the arrayIn by neighbering k using convolution,
% takes 1/3 of the time of AH_expandByK 
% if k=5, front expand 5, back expand k=5
% but has to be symmetrical expansion

if nargin < 2
    k = 15; %=0.0005*30000
end
if nargin < 3
    bound = [1,max(arrayIn)+k];
end
if isa(arrayIn, 'gpuArray') % if gpuArray, arrayAll has to be gpuArray too
    arrayAll = gpuArray(zeros(numel(bound),1));
else
    arrayAll = zeros(diff(bound)+1,1);
end
arrayAll(arrayIn) = 1;
mask = ones(2*k+1,1);
arrayConv = conv(arrayAll,mask,'same');
arrayOut = find(arrayConv>=1);
if nargin >= 3 % if bound exists
arrayOut(arrayOut<bound(1) | arrayOut>bound(2)) = []; % delete elements outside bound
end