function pxx = ns_fftpower(ts)
% function calculates power with a hann window

% apply a window
w = window(@hann,length(ts));
ts_for_fft = bsxfun(@times, ts, w);

% fft power
pxx = abs(fft(ts_for_fft)).^2;