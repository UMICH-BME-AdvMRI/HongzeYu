% HW3 Pr2 BIOMEDE 599-020

% Part a
data = load('Data_Assignment3_Problem2.mat');
coil_maps = data.coilmaps;
ks = data.kspaceData;
[npe, nro, ncoil] = size(ks);

image_full = zeros(npe, nro, ncoil);
coil_combined_image = zeros(npe, nro);
figure;
for idx=1:ncoil
    image_full(:,:,idx) = ifftshift(ifft2(ks(:,:,idx)));
    coil_combined_image = coil_combined_image + conj(coil_maps(:,:,idx)) .* image_full(:,:,idx);

    subplot(4,4,idx), imshow(abs(image_full(:,:,idx)), []), title(['Coil image ' num2str(idx)]);
    subplot(4,4,idx+ncoil), imshow(abs(coil_maps(:,:,idx)), []), title(['Coil map ' num2str(idx)]);
end
saveas(gcf, 'figure/pr2/part_a_coil_image_map.png')

figure; 
imshow(abs(coil_combined_image), []);
title('Coil combined image magnitude');
saveas(gcf, 'figure/pr2/part_a_coil_combined.png')

% part b    R=2
ks_R2 = zeros(npe, nro, ncoil);
ks_R2(1:2:end,:,:) = ks(1:2:end,:,:);

image_R2 = zeros(npe, nro, ncoil);
coil_combined_image_R2 = zeros(npe, nro);
figure;
for idx=1:ncoil
    image_R2(:,:,idx) = ifftshift(ifft2(ks_R2(:,:,idx)));
    coil_combined_image_R2 = coil_combined_image_R2 + conj(coil_maps(:,:,idx)) .* image_R2(:,:,idx);

    subplot(4,4,idx), imshow(abs(image_R2(:,:,idx)), []), title(['Coil image ' num2str(idx)]);
    subplot(4,4,idx+ncoil), imshow(abs(coil_maps(:,:,idx)), []), title(['Coil map ' num2str(idx)]);
end
saveas(gcf, 'figure/pr2/part_b_R2_coil_image_map.png')

figure; 
imshow(abs(coil_combined_image_R2), []);
title('Coil combined image magnitude R=2');
saveas(gcf, 'figure/pr2/part_b_R2.png')

% part c    SENSE R=2 
image_SENSE_R2 = zeros(npe, nro);

% solve y = Sx for x => x = S^-1 y
for ro=1:nro
    for pe=1:npe/2
        y = image_R2(pe, ro,:);
        y = reshape(y, [ncoil,1]); % aliased pixel across coil
        S = [];
        for coil=1:ncoil
            S = [S; coil_maps(pe, ro, coil) coil_maps(pe + npe/2, ro, coil)]; % coil sensitivity matrix
        end
        x = pinv(S) * y;
        image_SENSE_R2(pe,ro) = x(1);
        image_SENSE_R2(pe+npe/2, ro) = x(2);
    end
end

diff_image_R2 = coil_combined_image - 2 * image_SENSE_R2;
figure;
subplot(1,2,1), imshow(abs(image_SENSE_R2), []), title("SENSE recon R=2");
subplot(1,2,2), imshow(abs(diff_image_R2), []), title("diff image R=2");
saveas(gcf, 'figure/pr2/part_c_SENSE_image_and_diff_R2.png')

% part d    SENSE R=4
ks_R4 = zeros(npe, nro, ncoil);
ks_R4(1:4:end,:,:) = ks(1:4:end,:,:);
image_R4 = zeros(npe, nro, ncoil);

figure;
for idx=1:ncoil
    image_R4(:,:,idx) = ifftshift(ifft2(ks_R4(:,:,idx)));
    subplot(4,4,idx), imshow(abs(image_R4(:,:,idx)), []), title(['Coil image ' num2str(idx)]);
    subplot(4,4,idx+ncoil), imshow(abs(coil_maps(:,:,idx)), []), title(['Coil map ' num2str(idx)]);
end
saveas(gcf, 'figure/pr2/part_d_coil_image_map_R4.png')

image_SENSE_R4 = zeros(npe, nro);

% solve y = Sx for x => x = S^-1 y
for ro=1:nro
    for pe=1:npe/4
        y = image_R4(pe, ro,:);
        y = reshape(y, [ncoil,1]); % aliased pixel across coil
        S = [];
        for coil=1:ncoil
            S = [S; coil_maps(pe, ro, coil) coil_maps(pe + npe/4, ro, coil) coil_maps(pe + npe/2, ro, coil) coil_maps(pe + npe*3/4, ro, coil)]; % coil sensitivity matrix
        end
        x = pinv(S) * y;
        image_SENSE_R4(pe,ro) = x(1);
        image_SENSE_R4(pe+npe/4, ro) = x(2);
        image_SENSE_R4(pe+npe/2, ro) = x(3);
        image_SENSE_R4(pe+npe*3/4, ro) = x(4);
    end
end

diff_image_R4 = coil_combined_image - 4 * image_SENSE_R4;
figure;
subplot(1,2,1), imshow(abs(image_SENSE_R4), []), title("SENSE recon R=4");
subplot(1,2,2), imshow(abs(diff_image_R4), []), title("diff image R=4");
saveas(gcf, 'figure/pr2/part_d_SENSE_image_and_diff_R4.png')

close all;
