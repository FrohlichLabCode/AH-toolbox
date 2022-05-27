function mat2 = AH_resample(mat1,newFs,oldFs) % mat is freq x time

mat2 = resample(mat1',newFs,oldFs)'; % resample treats each column of x as an independent channel
mat2(:,1) = mat2(:,2); mat2(:,size(mat2,2)) = mat2(:,size(mat2,2)-1); % edge values are weird so just change to their neighbor value
% smooth out the edge effect, tried padding before resample, doesn't work
% as good, so just replace the 2 edge values with their neighbors
end