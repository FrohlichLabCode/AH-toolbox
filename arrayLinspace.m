function array = arrayLinspace(a,b,n,optionallinOrlog)
if nargin > 3
    linOrlog = optionallinOrlog;
else
    linOrlog = 1;
end
% a, b must be 1D array of the same size
a = reshape(squeeze(a),1,[]);
b = reshape(squeeze(b),1,[]);
array = NaN(n,numel(a));
for i = 1:numel(a)
    if linOrlog == 0
        array(:,i) = linspace(a(i),b(i),n);
    else
        array(:,i) = exp(linspace(log(a(i)),log(b(i)),n));
    end
end
return
end