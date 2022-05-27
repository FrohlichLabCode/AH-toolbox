function spec = nan2one(spec)
% The function takes into a spectrogram and convert NaN and Inf into 1 (so
% that when take pow2db(1)=0.
% Angel 2020/4/23

mask = isnan(spec) | isinf(spec);
spec(mask) = 1;
end