% HW2 Pr2

%   load the maps
tissue_maps = load('brain_maps.mat');
T1map = tissue_maps.T1map;
T2map = tissue_maps.T2map;
M0map = tissue_maps.M0map;
image_size = size(T1map);

%% Pr2(a) single-echo SE

contrast = 'PD'; % 'T1', 'T2' or 'PD'

switch contrast
    case 'T1'
        %   T1-weighted image
        %   For T1w images, choose short TE and TR
        TE = 15; % ms
        TR = 500; % ms
    case 'T2'
        %   T2-weighted image
        %   For T2w images, choose long TE and TR
        TE = 100; % ms
        TR = 3000; % ms 
    case 'PD'
        %   PD-weighted image
        %   For PDw images, choose short TE and long TR
        TE = 15; % ms
        TR = 5000; % ms 
end
folderPath = 'plot/pr2_a';

%   single-echo SE

%   Implementation:
%   As we only acqure one kspace line during one echo, we need to simulate 
%   in total number of kspace lines of echo. (same for FSE)
%   And for each echo, we use EPG to simulate the signal, take FT to get
%   its kspace, and then take the corresponding kspace line. Eventually,
%   combine total 256 lines from 256 echoes together to form the acquired 
%   kspace and use IFT to get the image. Thus in total 256 TR.
%   E.g. for the nth echo (n=1,...,256), we should take the nth ks line

acquired_ks = zeros(image_size);
% image series, 256 kspace lines mean 256 echo images
num_echoes = 256;
image = zeros(image_size(1), image_size(2), num_echoes);
n = 0;

for xx=1:image_size(1)
    for yy=1:image_size(2)
        % simulate each pixel one-by-one
        m = [0 0 1]';
        for echo_num=1:num_echoes
            [m, sig] = EPG_single_echo_SE(m,T1map(xx,yy),T2map(xx,yy),TE,TR);
            image(xx,yy,echo_num) = sig * M0map(xx,yy); % apply PD map for each sig
            % note PD is only applied for each pixel at each echo once
        end
        n = n + 1;
        if mod(n, 1000) == 0
            fprintf("The %i pixel finished\n", n);
        end
    end
end

% The code is able to reduce the time to implement an EPI-like acquisition
% by changing the num_echoes, by setting num_echoes = 256, we are acquire
% one line per echo by default as required in HW
for echo_num=1:num_echoes
    echo_image = squeeze(image(:,:,echo_num));
    echo_image(isnan(echo_image)) = 0;
    ks = fftshift(fft2(echo_image));

    % for the sake of running time, assume Nlines are acquired during one
    % echo like EPI
    Nline = image_size(1) / num_echoes;
    start_line = (echo_num-1) * Nline + 1;
    end_line = echo_num * Nline;
    acquired_ks(start_line:end_line, :) = ks(start_line:end_line, :);
end

acquired_image = ifft2(ifftshift(acquired_ks));
imshow(abs(acquired_image), []);
title([contrast 'w image']);
saveas(gcf, fullfile(folderPath, [contrast 'w.png']));




%% Pr2(b) FSE

TR = 3000; % ms     TR = 3s
ESP = 5; % ms
ETL = 32;

% part i
tissue_pair = {[1000, 50], [1000, 100], [2000, 50], [2000, 100]};
tt = 0:0.5:TR;
figure;
legends = cell(1, length(tissue_pair));

for i = 1:length(tissue_pair)
    transverse_m = [];
    T1 = tissue_pair{i}(1);
    T2 = tissue_pair{i}(2);
    m = [0 0 1]';

    for idx=1:4 % first 4 TR
        [m, sig] = EPG_FSE(m, T1, T2, TR, ESP, ETL);
    end

    % Last TR - dumb implementation
    m = EPG_RF(m, pi/2, 0); % 90 x pulse
    transverse_m = [transverse_m abs(m(1,1))]; % t = 0

    for echo_num = 1:ETL
        for num = 1:4
            m = EPG_relax(m, T1, T2, 0.5);
            transverse_m = [transverse_m abs(m(1,1))]; % t = 0.5 - 2.0
        end

        m = EPG_grad(m,1);
        m = EPG_RF(m, pi, pi/2);
        m = EPG_grad(m,1);

        transverse_m = [transverse_m abs(m(1,1))]; % t = 2.5

        for num = 1:5
            m = EPG_relax(m, T1, T2, 0.5);
            transverse_m = [transverse_m abs(m(1,1))]; % t = 3.0 - 5.0
        end
    end
    % t = 180
    num_samples_left = TR/0.5 - ETL * ESP/0.5;
    for num=1:num_samples_left
        m = EPG_relax(m, T1, T2, 0.5);
        transverse_m = [transverse_m abs(m(1,1))];
    end

    plot(tt, transverse_m);
    hold on;
    legends{i} = sprintf('T1 = %d, T2 = %d', T1, T2);
end

legend(legends);
title("Transverse magnetization last TR")
xlabel('time (ms)')
ylabel('Trans magnetization')

savepath = 'plot/pr2_b_i';
saveas(gcf, fullfile(savepath, 'Transverse_Magnetization.png'));


% part ii
num_of_TR = image_size(1) / ETL; % How many TR needed (how many group of echo train

% prepare the kspace
acquired_ks = zeros(image_size(1), image_size(2));

% for each TR, there are ETL number of magnetization
echo_image = zeros(image_size(1), image_size(2), image_size(1));

n = 0;
% simulate 256 echo image for later kspace acquisition
for xx=1:image_size(1)
    for yy=1:image_size(2)
        m = [0 0 1]';
        for TR_num = 1:num_of_TR % for each TR
            [m, sig] = EPG_FSE(m, T1map(xx,yy), T2map(xx,yy), TR, ESP, ETL);
            PD_weighting = M0map(xx,yy);
            % sig - size (1,ETL)
            sig = sig * PD_weighting; % apply PD weighting

            for idx=1:ETL
                echo_image(xx,yy,idx+ETL*(TR_num-1))=sig(idx);
            end
        end

        % count currently on which pixel
        n = n + 1;
        if mod(n, 1000) == 0
            fprintf("The %i pixel finished\n", n);
        end
    end
end

% Phase-encoding interleaving acquisition
% assume the very center of kspace is 256/2 = 128
% as here Teff is 80, the 16 line and 17 line should be in the ks center
Teff = 80; % ms
eff_echo_idx = Teff / ESP; % The echo idx for Teff

for TR_num = 1:num_of_TR
    % for example, for the first TR, if ETL = 32, we acquire
    % line 1,9, ...., 121, 129, ..., 246
    % the 16th echo is line 121, 17 line is 129
    lines_to_acquire = TR_num:num_of_TR:image_size(1);
    for echo_idx=1:length(lines_to_acquire)
        % which ks line to acquire
        current_line = lines_to_acquire(echo_idx);
        % which echo are we currently at
        current_echo = (TR_num-1)*ETL + echo_idx;
        echo_magnetization = echo_image(:,:, current_echo);
        echo_magnetization(isnan(echo_magnetization)) = 0;

        echo_ks = fftshift(fft2(echo_magnetization));
        acquired_ks(current_line, :) = echo_ks(current_line,:);
    end
end

acquired_image = ifft2(ifftshift(acquired_ks));
imshow(abs(acquired_image));
titleStr = sprintf("ETL = %d, Teff= %d", ETL, Teff);
title(titleStr);
folderPath = 'plot/pr2_b_ii';
filename = sprintf("ETL_%d_Teff_%d.png", ETL, Teff);
saveas(gcf, fullfile(folderPath, filename));

% showing the kspace fill order
kspace_diagram = zeros(image_size(1), image_size(2));

for TR_num=1:num_of_TR
    kspace_diagram(TR_num:num_of_TR:end, :) = TR_num;
end

figure;
imshow(mat2gray(kspace_diagram))
title('K-space Sampling Order');
xlabel('Readout');
ylabel('Phase Encoding');
axis on;
grid on; 
folderPath = 'plot/pr2_b_ii';
filename = sprintf("ks_sampling_order.png");
saveas(gcf, fullfile(folderPath, filename))

% part iii
% echo image is the same, we don't need to do extra EPG simulation here
% as ETL, ESP not changed
% only need to rewrite the k-space acquisition part
Teff = 120; % 40 or 120
echo_idx = Teff / ESP; % The 8th or 24th echo is kspace center
shift_idices = Teff / ESP - ETL / 2;

% Phase-encoding interleaving acquisition
% assume the very center of kspace is 256/2 = 128
% as here Teff is 40, the 8th line and 9th line should be in the ks center

for TR_num = 1:num_of_TR
    % for example, for the first TR, if ETL = 32, we acquire
    % line 1,9, ...., 121, 129, ..., 246
    % the 16th echo is line 121, 17 line is 129
    lines_to_acquire = TR_num:num_of_TR:image_size(1);
    lines_to_acquire = circshift(lines_to_acquire, shift_idices); % make the 8th echo at 16th line
    for echo_idx=1:length(lines_to_acquire)
        % which ks line to acquire
        current_line = lines_to_acquire(echo_idx);
        % which echo are we currently at
        current_echo = (TR_num-1)*ETL + echo_idx;
        echo_magnetization = echo_image(:,:, current_echo);
        echo_magnetization(isnan(echo_magnetization)) = 0;

        echo_ks = fftshift(fft2(echo_magnetization));
        acquired_ks(current_line, :) = echo_ks(current_line,:);
    end
end

acquired_image = ifft2(ifftshift(acquired_ks));
imshow(abs(acquired_image));
titleStr = sprintf("ETL = %d, Teff= %d", ETL, Teff);
title(titleStr);
folderPath = 'plot/pr2_b_iii';
filename = sprintf("ETL_%d_Teff_%d.png", ETL, Teff);
saveas(gcf, fullfile(folderPath, filename));

% kspace filling order plot are using the same plot of previos question.

% part iv

Teff = 80;
TR = 3000; % ms     TR = 3s
ESP = 5; % ms
ETL = 64; % 16, 64, 128
eff_echo_idx = Teff / ESP; % The echo idx for Teff = 16
shift_idices = eff_echo_idx - ETL / 2; % 8 for 16, -16 for 64, -48 for 128 


num_of_TR = image_size(1) / ETL; % How many TR needed (how many group of echo train

% prepare the kspace
acquired_ks = zeros(image_size(1), image_size(2));

% for each TR, there are ETL number of magnetization
echo_image = zeros(image_size(1), image_size(2), image_size(1));

n = 0;
% simulate 256 echo image for later kspace acquisition
for xx=1:image_size(1)
    for yy=1:image_size(2)
        m = [0 0 1]';
        for TR_num = 1:num_of_TR % for each TR
            [m, sig] = EPG_FSE(m, T1map(xx,yy), T2map(xx,yy), TR, ESP, ETL);
            PD_weighting = M0map(xx,yy);
            % sig - size (1,ETL)
            sig = sig * PD_weighting; % apply PD weighting

            for idx=1:ETL
                echo_image(xx,yy,idx+ETL*(TR_num-1))=sig(idx);
            end
        end

        % count currently on which pixel
        n = n + 1;
        if mod(n, 1000) == 0
            fprintf("The %i pixel finished\n", n);
        end
    end
end

% Phase-encoding interleaving acquisition
% assume the very center of kspace is 256/2 = 128
% as here Teff is 80, the 16 line should be in the ks center

for TR_num = 1:num_of_TR
    % for example, for the first TR, if ETL = 32, we acquire
    % line 1,9, ...., 121, 129, ..., 246
    % the 16th echo is line 121, 17 line is 129
    lines_to_acquire = TR_num:num_of_TR:image_size(1);
    lines_to_acquire = circshift(lines_to_acquire, shift_idices);
    for echo_idx=1:length(lines_to_acquire)
        % which ks line to acquire
        current_line = lines_to_acquire(echo_idx); % echo idx within this TR
        % which echo are we currently at
        current_echo = (TR_num-1)*ETL + echo_idx; % total idx for this echo
        echo_magnetization = echo_image(:,:, current_echo);
        echo_magnetization(isnan(echo_magnetization)) = 0;

        echo_ks = fftshift(fft2(echo_magnetization));
        acquired_ks(current_line, :) = echo_ks(current_line,:);
    end
end

acquired_image = ifft2(ifftshift(acquired_ks));
imshow(abs(acquired_image));
titleStr = sprintf("ETL = %d, Teff= %d", ETL, Teff);
title(titleStr);
folderPath = 'plot/pr2_b_iv';
filename = sprintf("ETL_%d_Teff_%d.png", ETL, Teff);
saveas(gcf, fullfile(folderPath, filename));


%   single-echo SE

%   Implementation:
%   As we only acqure one kspace line during one echo, we need to simulate 
%   in total number of kspace lines of echo. (same for FSE)
%   And for each echo, we use EPG to simulate the signal, take FT to get
%   its kspace, and then take the corresponding kspace line. Eventually,
%   combine total 256 lines from 256 echoes together to form the acquired 
%   kspace and use IFT to get the image. Thus in total 256 TR.
%   E.g. for the nth echo (n=1,...,256), we should take the nth ks line
TE = 80;
TR = 3000;
acquired_ks = zeros(image_size);
% image series, 256 kspace lines mean 256 echo images
num_echoes = 256;
image = zeros(image_size(1), image_size(2), num_echoes);
n = 0;

for xx=1:image_size(1)
    for yy=1:image_size(2)
        % simulate each pixel one-by-one
        m = [0 0 1]';
        for echo_num=1:num_echoes
            [m, sig] = EPG_single_echo_SE(m,T1map(xx,yy),T2map(xx,yy),TE,TR);
            image(xx,yy,echo_num) = sig * M0map(xx,yy); % apply PD map for each sig
            % note PD is only applied for each pixel at each echo once
        end
        n = n + 1;
        if mod(n, 1000) == 0
            fprintf("The %i pixel finished\n", n);
        end
    end
end

% The code is able to reduce the time to implement an EPI-like acquisition
% by changing the num_echoes, by setting num_echoes = 256, we are acquire
% one line per echo by default as required in HW
for echo_num=1:num_echoes
    echo_image = squeeze(image(:,:,echo_num));
    echo_image(isnan(echo_image)) = 0;
    ks = fftshift(fft2(echo_image));

    % for the sake of running time, assume Nlines are acquired during one
    % echo like EPI
    Nline = image_size(1) / num_echoes;
    start_line = (echo_num-1) * Nline + 1;
    end_line = echo_num * Nline;
    acquired_ks(start_line:end_line, :) = ks(start_line:end_line, :);
end

folderPath = 'plot/pr2_b_iv';
acquired_image = ifft2(ifftshift(acquired_ks));
imshow(abs(acquired_image), []);
title("Single echo SE TE = 80 ms");
saveas(gcf, fullfile(folderPath, 'single_echo_SE.png'));



















