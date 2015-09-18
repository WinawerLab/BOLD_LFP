function alpha_signal = ns_alpha_signal(alpha_inputs,poisson_rate_a,dt,do_plot)

srate = 1/dt;

%%%% DESIGN ALPHA FILTER
band = [8 15];
a_Rp = 3; a_Rs = 60; % order Butterworth
a_delta = 0.001*2/srate;
a_low_p = band(2)*2/srate;
a_high_p = band(1)*2/srate;
a_high_s = max(a_delta, a_high_p - 0.1);
a_low_s = min(1-a_delta, a_low_p + 0.1);
[a_n_band, a_wn_band] = buttord([a_high_p a_low_p], [a_high_s a_low_s], a_Rp, a_Rs);
[a_bf_b, a_bf_a] = butter(a_n_band, a_wn_band);

%%%% DESIGN LOWPASS FILTER
band = [3];
low_Rp = 3; low_Rs = 60; % order Butterworth
low_high_p = band(1)*2/srate;
low_high_s = (band(1)+20)*2/srate;
[low_n_band, low_wn_band] = buttord(low_high_p, low_high_s, low_Rp, low_Rs);
[low_bf_b, low_bf_a] = butter(low_n_band, low_wn_band,'low');

alpha_signal = zeros(size(alpha_inputs));

% loop all alpha neurons
for k=1:size(alpha_inputs,2)
    % get alpha signal to use:
    alpha_use = alpha_inputs(:,k);

    % filter for alpha
%     alpha_useF = poisson_rate_a*filtfilt(a_bf_b, a_bf_a, alpha_use);
    
%     %%%% option 1, but alpha does not correlate with mean response
%     %%%% anymore...
%     % add (DC offset) to signal
%     alpha_signal(:,k) = alpha_useF - 1/(1+poisson_rate_a);
%     
%     % make sure there is no broadband added with adding the envelope
%     alpha_useF_fft = fft(alpha_useF);
%     alpha_signal_fft = fft(alpha_signal(:,k));
%     alpha_signal_fft(10:end-9) = alpha_useF_fft(10:end-9);
%     alpha_signal(:,k) = real(ifft(alpha_signal_fft));

    %%%% option 2, complex...
%     % zero-pad alpha inputs
%     alpha_use = cat(1,zeros(size(alpha_use)),alpha_use,zeros(size(alpha_use)));
% 
%     % filter for alpha
%     alpha_useF = poisson_rate_a*filtfilt(a_bf_b, a_bf_a, alpha_use);
%     
%     % get alpha envelope
%     alpha_envelope = abs(hilbert(alpha_useF));
%     
%     % invert envelope
%     alpha_envelope = -1./(1+alpha_envelope);
% 
%     % filter the envelope
%     alpha_envelopeF = filtfilt(low_bf_b, low_bf_a, alpha_envelope);
%     
%     % get non-zeropadded alpha signal and envelope back
%     alpha_useFS = alpha_useF(size(alpha_inputs,1)+1:size(alpha_inputs,1)*2);
%     alpha_envelopeFS = alpha_envelopeF(size(alpha_inputs,1)+1:size(alpha_inputs,1)*2);
%            
%     % add filtered envelope (DC offset) to signal
%     alpha_signal(:,k) = alpha_useFS + alpha_envelopeFS;
%     
%     % make sure there is no broadband added with adding the envelope
%     alpha_useFS_fft = fft(alpha_useFS);
%     alpha_signal_fft = fft(alpha_signal(:,k));
%     alpha_signal_fft(10:end-9) = alpha_useFS_fft(10:end-9);
%     alpha_signal(:,k) = real(ifft(alpha_signal_fft));
    
    %%%% option 3, 
    % add assymetry to signal
    % zero-pad alpha inputs
    alpha_use = cat(1,zeros(size(alpha_use)),alpha_use,zeros(size(alpha_use)));

    % filter for alpha
    alpha_useF = poisson_rate_a*filtfilt(a_bf_b, a_bf_a, alpha_use);
    
    % get alpha envelope
    alpha_envelope = abs(hilbert(alpha_useF));
    
    % filter the envelope
    alpha_envelopeF = filtfilt(low_bf_b, low_bf_a, alpha_envelope);
    
    % get non-zeropadded alpha signal and envelope back
    alpha_useFS = alpha_useF(size(alpha_inputs,1)+1:size(alpha_inputs,1)*2);
    alpha_envelopeFS = alpha_envelopeF(size(alpha_inputs,1)+1:size(alpha_inputs,1)*2);
           
    % add filtered envelope (DC offset) to signal
    alpha_signal(:,k) = alpha_envelopeFS.*(1+alpha_useFS);

    % make sure there is no broadband added with adding the envelope
    alpha_useFS_fft = fft(alpha_useFS);
    alpha_signal_fft = fft(alpha_signal(:,k));
    alpha_signal_fft(10:end-9) = alpha_useFS_fft(10:end-9);
    alpha_signal(:,k) = real(ifft(alpha_signal_fft));

    %%%%% plot if do_plot==1 for k=1
    if k==1 && do_plot==1
        figure
        subplot(3,2,[1:2]),hold on
        plot(alpha_useF,'k')
        xlim([0 length(alpha_use)])
        plot(alpha_envelope,'c')
        plot(alpha_envelopeF,'b--')

        subplot(3,2,[3 5]),hold on
        plot(alpha_useFS,'k')
        plot(alpha_envelopeFS,'c')
        plot(alpha_signal(:,k),'r') 
        
        subplot(3,2,[4 6]),hold on
        [pxx,f] = pwelch(alpha_useFS,srate,0,srate,srate);
        plot(f,log10(pxx),'k')
        [pxx,f] = pwelch(alpha_envelopeFS,srate,0,srate,srate);
        plot(f,log10(pxx),'c')
        [pxx,f] = pwelch(alpha_signal(:,k),srate,0,srate,srate);
        plot(f,log10(pxx),'r')
        xlim([0 200])
    end
    clear alpha_use alpha_envelope
end
