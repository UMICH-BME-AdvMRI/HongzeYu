% HW3 Pr1 BIOMEDE 599-020

% Part a
kspace_data = load("Data_Assignment3_Problem1.mat");
ks = kspace_data.kspaceData_SingleCoil;

% Retrospective undersampling
[npe, nro] = size(ks);
ks_zerofilled = zeros(npe, nro);
ks_zerofilled(1:(5 * npe/8), :) = ks(1:(5 * npe/8), :);

% Partial Fourier Image
image_undersampled = ifftshift(ifft2(ks_zerofilled));
mag_undersampled = abs(image_undersampled);
phase_undersampled = angle(image_undersampled);
figure;
subplot(1, 2, 1), imshow(mag_undersampled, []), title('Magnitude - Partial Fourier');
subplot(1, 2, 2), imshow(phase_undersampled, []), title('Phase - Partial Fourier');
saveas(gcf, 'figure/pr1/part_a_partial_Fourier.png')

% For comparison
image_full = ifftshift(ifft2(ks));
diff_image = image_full - image_undersampled;
mag_diff = abs(diff_image);
phase_diff = angle(diff_image);
figure;
subplot(1, 2, 1), imshow(mag_diff, []), title('Magnitude Difference');
subplot(1, 2, 2), imshow(phase_diff, []), title('Phase Difference');
saveas(gcf, 'figure/pr1/part_a_partial_Fourier_diff.png')

% part b - POCS

% Phase estimation
ks_lowres = zeros(npe, nro);
ks_lowres(3*npe/8:5*npe/8, :) = ks(3*npe/8:5*npe/8, :);

% define hanning filter
hanning_filter = hann(npe) * ones(1, nro);
ks_lowres = ks_lowres .* hanning_filter;
image_lowres = ifftshift(ifft2(ks_lowres));
phase_lowres = angle(image_lowres);

% iteative reconstruction
ks_iterative = ks_zerofilled; % initialization 
for idx=1:10
    image_iterative = ifftshift(ifft2(ks_iterative));
    image_iterative = abs(image_iterative) .* exp(1i * phase_lowres);
    ks_iterative = fft2(fftshift(image_iterative));
    ks_iterative(1:(5 * npe/8), :) = ks_zerofilled(1:(5 * npe/8), :);
end

% Display POCS reconstruction
image_pocs = ifftshift(ifft2(ks_iterative));
mag_pocs = abs(image_pocs);
phase_pocs = angle(image_pocs);
figure;
subplot(1, 2, 1), imshow(mag_pocs, []), title('Magnitude - POCS');
subplot(1, 2, 2), imshow(phase_pocs, []), title('Phase - POCS');
saveas(gcf, 'figure/pr1/part_b_POCS.png')

diff_pocs_image = image_full - image_pocs;
mag_diff_pocs = abs(diff_pocs_image);
phase_diff_pocs = angle(diff_pocs_image);
figure;
subplot(1, 2, 1), imshow(mag_diff_pocs, []), title('POCS Magnitude Difference');
subplot(1, 2, 2), imshow(phase_diff_pocs, []), title('POCS Phase Difference');
saveas(gcf, 'figure/pr1/part_b_POCS_diff.png')

close all;

