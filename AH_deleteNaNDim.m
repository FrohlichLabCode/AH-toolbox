function deleteMask = AH_getNaNDimMask(mat,otherDim)
% delete NaN or 0 sessions
nanMask = all(isnan(mat),otherDim);
zeroMask = all(mat==0,otherDim);
deleteMask = nanMask | zeroMask;

                