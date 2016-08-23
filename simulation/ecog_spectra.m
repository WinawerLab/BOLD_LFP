function [pxx,f] = ecog_spectra(ts,srate,nr_freq)
% function calculates power with a hann window
% nr_freq = number of frequencies to return

w = window(@hann,srate/4);
[pxx,f] = pwelch(ts,w,length(w)/2,srate,srate);

f = f(1:nr_freq);
pxx = pxx(1:nr_freq);
