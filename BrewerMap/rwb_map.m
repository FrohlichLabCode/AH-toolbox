function c = rwb_map
% This shorthand script calls for brewermap and get its red-white-blue map
c = flipud(brewermap([],'RdBu')); % default 64 gradients
% c = flipud(brewermap([],'Spectral')) % this gives red-yellow-blue map
end