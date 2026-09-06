% Example usage of the blind_deconvolution function

% Generate an example signal and PSF
original_signal = sin(2 * pi * (0:0.01:1));
psf = exp(-(-10:10).^2 / 2);

% Convolve the signal with the PSF to simulate the observed signal
observed_signal = conv(original_signal, psf, 'same');

% Initial guess for the PSF (random initialization)
initial_psf = randn(size(psf));

% Perform blind deconvolution
max_iter = 100;
[recovered_signal, recovered_psf] = blind_deconvolution(observed_signal, initial_psf, max_iter);

% Plot the results
figure;
subplot(3, 1, 1);
plot(original_signal);
title('Original Signal');

subplot(3, 1, 2);
plot(observed_signal);
title('Observed Signal');

subplot(3, 1, 3);
plot(recovered_signal);
title('Recovered Signal');
