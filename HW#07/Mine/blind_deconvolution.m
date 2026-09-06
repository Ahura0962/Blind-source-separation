function [recovered_signal, recovered_psf] = blind_deconvolution(observed_signal, initial_psf, max_iter)
    % Function to perform blind deconvolution in the frequency domain
    % observed_signal: The observed (blurred) signal
    % initial_psf: Initial estimate of the point spread function (PSF)
    % max_iter: Maximum number of iterations

    % Convert the observed signal and initial PSF to frequency domain
    observed_signal_fft = fft(observed_signal);
    psf_fft = fft(initial_psf, length(observed_signal)); % Pad PSF to the same length

    % Initialize the recovered signal
    recovered_signal_fft = observed_signal_fft ./ psf_fft;

    for iter = 1:max_iter
        % Update the estimated signal
        recovered_signal = real(ifft(recovered_signal_fft));
        
        % Recalculate the PSF in the frequency domain
        estimated_signal_fft = fft(recovered_signal);
        psf_fft = observed_signal_fft ./ estimated_signal_fft;
        
        % Update the recovered signal in the frequency domain
        recovered_signal_fft = observed_signal_fft ./ psf_fft;
        
        % Ensure the PSF is real and non-negative
        recovered_psf = real(ifft(psf_fft));
        recovered_psf(recovered_psf < 0) = 0;
        
        % Update the PSF in the frequency domain
        psf_fft = fft(recovered_psf, length(observed_signal));
    end
    
    % Convert the final results back to time domain
    recovered_signal = real(ifft(recovered_signal_fft));
    recovered_psf = real(ifft(psf_fft));
end
