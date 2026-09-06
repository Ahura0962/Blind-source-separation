function [signal_normalized_f, norm_sig] = Frequency_Normalization(signal_f)
    signal_t = ifft(ifftshift(signal_f));
    norm_sig = vecnorm(signal_t);
    signal_t = signal_t ./ norm_sig;
    signal_normalized_f = fftshift(fft(signal_t));
end

