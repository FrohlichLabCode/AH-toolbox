% compute hilbert transform on lfp input; Fc are the filter bounds
% if have predefined bounds, input Freq as [lowfreq highfreq] ; otherwise,
% input as center frequency
function C = cz_computeHilbert(lfpInput, Freq, Fs)

if numel(Freq) == 1 % if bounds are already defined
    hilBounds = 0.1; % percent above below center freq for hilbert filtr; 0.05 gives minimal overlap
    Fc = [Freq-hilBounds*Freq Freq+hilBounds*Freq]; % Cut-off frequencies
else
    Fc = Freq;
end

Nyq = Fs/2;
Wn = Fc/Nyq; % divide by nyquist
[B, A] = butter(2, Wn); % calc filter parameters

C_tmp = hilbert(filtfilt(B, A, lfpInput')); %filter and hilbert; transpose b/c functions apply to columns
C = C_tmp';

end