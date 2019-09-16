%An improved method for R-peak detection by using Shannon energy envelope
%An R-peak detection method based on peaks of Shannon energy envelope
function [out_sig,peak_locs_final] = ISEE(signal,fs)
    % input: original signal, sample frequency
    % output: filtered signal, detected peak locations
    % bandpass filter
    BP_filter = chebyshevI_bandpass(4,fs,6,18,'cheby1');
    % forward filtering
    f_signal_forward = filter(BP_filter,signal);
    % backward filtering
    f_signal = fliplr(f_signal_forward);
    f_signal = filter(BP_filter,f_signal);
    f_signal = fliplr(f_signal);
    %% this part checks whether the filter is working
%     sig_mags = abs(fft(signal));
%     f_sig_mags = abs(fft(f_signal));
%     f_sigf_mags = abs(fft(f_signal_forward));
%     N = length(signal);
%     
%     figure,
%     subplot(3,1,1)
%     bin_vals = [0 : N-1];
%     fax_Hz = bin_vals*fs/N;
%     N_2 = ceil(N/2);
%     plot(fax_Hz(2:N_2), sig_mags(2:N_2))
%     
%     subplot(3,1,2)
%     bin_vals = [0 : N-1];
%     fax_Hz = bin_vals*fs/N;
%     N_2 = ceil(N/2);
%     plot(fax_Hz(2:N_2), f_sigf_mags(2:N_2))
%     
%     subplot(3,1,3)
%     bin_vals = [0 : N-1];
%     fax_Hz = bin_vals*fs/N;
%     N_2 = ceil(N/2);
%     plot(fax_Hz(2:N_2), f_sig_mags(2:N_2))
    %% other filters for the detection
    % first order differentiation
    d_n = f_signal(2:end)-f_signal(1:end-1);
    % normalize the signal
    norm_dn = abs(d_n)/max(abs(d_n));
    % shannon energy envelope
    se_n = (-1)*(norm_dn.^2).*log(norm_dn.^2);
    % moving average filter
    N = 65;
    ma_filter = (1/N)*ones(1,N);
    s_n = filter(ma_filter,1,se_n);
%     s_n = conv(se_n,ma_filter,'same'); % another way of doing average
%     filtering
    % first order differentiation
    d_n_2 = s_n(2:end)-s_n(1:end-1);
    % normalize the signal
    norm_sig_2 = abs(d_n_2)/max(abs(d_n_2));
    % squaring
    s_n_2 = norm_sig_2.^2;
    % moving average filter
    N = 85;
    ma_filter = (1/N)*ones(1,N);
    out_sig = filter(ma_filter,1,s_n_2);

    %% plot to see the output of each step    
    %     figure,
    %     subplot(10,1,1)
    %     plot(signal)
    %     subplot(10,1,2)
    %     plot(f_signal)
    %     subplot(10,1,3)
    %     plot(d_n)
    %     subplot(10,1,4)
    %     plot(norm_dn)
    %     subplot(10,1,5)
    %     plot(se_n)
    %     subplot(10,1,6)
    %     plot(s_n)
    %     subplot(10,1,7)
    %     plot(d_n_2)
    %     subplot(10,1,8)
    %     plot(norm_sig_2)
    %     subplot(10,1,9)
    %     plot(s_n_2)
    %     subplot(10,1,10)
    %     plot(out_sig)

    %% find the peaks
    [~,peak_locs_temp]=findpeaks(out_sig);
    % detect real r peaks
    peak_locs_final = real_r_peak_detection_loop(f_signal,fs,peak_locs_temp,25);