%% Q3 Slice Profile Simulation

%% Q3.2

% Constants
gamma = 42.58*1e6;  % Hz/T
TBW = 8;        
fa = pi/2;  
slice_thickness = 5; % mm
tp = 2e-3;

% Time vector
t = [-1:0.001:1]*1e-3; % s 
dt = t(2)-t(1);

% Sinc RF pulse
BW = TBW / tp;  % Bandwidth of the pulse
rf_normalized = sinc(BW*t);
BW = BW * 2;

% Compute the necessary amplitude to achieve the desired flip angle
rf_area = sum(rf_normalized)*dt;
rf_amplitude = fa/(2*pi*gamma*rf_area)*100;  % Scaling factor for the sinc
rf = rf_amplitude * rf_normalized;

% Zero out the values after the pulse duration
rf = [rf, 0*rf]; 
rf_90 = rf;


tt = linspace(0, 4, length(rf));

% Plot RF
figure;
plot(tt, rf*1e4);
xlabel('Time (ms)');
ylabel('RF Amplitude (μT)');
title('RF waveform FA=90');
grid on;

% Calculate Gss
Gss = 18.8*1e-3; % T/m;

% Calculate ramp up and ramp down time for the gradient
Gz = Gss * ones(1, length(t));
Gz = [Gz, -Gz/2];


% % Plot Gz
figure;
plot(tt, Gz*1e3);
xlabel('Time (ms)');
ylabel('Gradient Amplitude mT/m');
title('Gradient Gz');
grid on;

% Constants for Bloch Simulation
T1 = 1000;  % ms
T2 = 100;   % ms option T2 = 2
% T2 = 2;
df = [0, 200]; % Off-resonance frequencies in Hz

% Compute the slice profile for each off-resonance frequency
position = [-1:0.01:1]*1e-2; % Position in mm
Msig_0 = zeros(1, length(position));
Msig_200 = zeros(1, length(position));
Mss_0 = zeros(3, length(position)); % 0Hz df
Mss_200 = zeros(3, length(position)); % 200Hz df

for freq_idx = 1:length(df)
    for j = 1: length(position)
        M = [0, 0, 1]';
        [A, B] = freeprecess(1000*dt/2, T1, T2, df(freq_idx));
        for k = 1:length(rf)
            M = A*M + B;
            phi = 2*pi*gamma*position(j)*Gz(k)*dt/2;
            M = zrot(phi)*M;
            M = throt(abs(rf(k)), angle(rf(k))) * M;
            M = A*M+B;
            M = zrot(phi)*M;
        end 

        if freq_idx == 1
            Mss_0(:, j) = M;
            Msig_0(j) = M(1) + 1i*M(2);
        end
        if freq_idx == 2
            Mss_200(:, j) = M;
            Msig_200(j) = M(1) + 1i*M(2);
        end
    end
end


figure
hold on; grid on;
plot(position*1e3, Mss_0(1, :))
plot(position*1e3, Mss_0(2, :))
plot(position*1e3, Mss_0(3, :))
plot(position*1e3, abs(Msig_0(:)))
legend('Mx', 'My', 'Mz', 'Mxy')
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization vs position off-resonace = 0Hz FA=90')

figure
hold on; grid on;
plot(position*1e3, Mss_200(1, :))
plot(position*1e3, Mss_200(2, :))
plot(position*1e3, Mss_200(3, :))
plot(position*1e3, abs(Msig_200(:)))
legend('Mx', 'My', 'Mz', 'Mxy')
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization vs position off-resonace = 200Hz FA=90')


%% Q3.3

% 30 degree
fa = 1/6 * pi;

rf_amplitude = fa/(2*pi*gamma*rf_area)*100;  % Scaling factor for the sinc
rf = rf_amplitude * rf_normalized;

% Zero out the values after the pulse duration
rf_30 = [rf, 0*rf]; 

% Compute the slice profile for each off-resonance frequency
position = [-1:0.01:1]*1e-2; % Position in mm
Msig_fa30 = zeros(1, length(position));
Mss_fa30 = zeros(3, length(position)); % 0Hz df
df = 0;


for j = 1: length(position)
    M = [0, 0, 1]';
    [A, B] = freeprecess(1000*dt/2, T1, T2, df);
    for k = 1:length(rf_30)
        M = A*M + B;
        phi = 2*pi*gamma*position(j)*Gz(k)*dt/2;
        M = zrot(phi)*M;
        M = throt(abs(rf_30(k)), angle(rf_30(k))) * M;
        M = A*M+B;
        M = zrot(phi)*M;
    end 
    Msig_fa30(j) = M(1) + 1i*M(2);
    Mss_fa30(:, j) = M;
end

figure
hold on; grid on;
plot(position*1e3, Mss_fa30(1, :))
plot(position*1e3, Mss_fa30(2, :))
plot(position*1e3, Mss_fa30(3, :))
plot(position*1e3, abs(Msig_fa30(:)))
legend('Mx', 'My', 'Mz', 'Mxy')
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization vs position FA=30 df=0')

% 10 degree
fa = 1/18 * pi;

rf_amplitude = fa/(2*pi*gamma*rf_area)*100;  % Scaling factor for the sinc
rf = rf_amplitude * rf_normalized;

% Zero out the values after the pulse duration
rf_10 = [rf, 0*rf]; 

% Compute the slice profile for each off-resonance frequency
position = [-1:0.01:1]*1e-2; % Position in mm
Msig_fa10 = zeros(1, length(position));
Mss_fa10 = zeros(3, length(position)); % 0Hz df
df = 0;


for j = 1: length(position)
    M = [0, 0, 1]';
    [A, B] = freeprecess(1000*dt/2, T1, T2, df);
    for k = 1:length(rf_10)
        M = A*M + B;
        phi = 2*pi*gamma*position(j)*Gz(k)*dt/2;
        M = zrot(phi)*M;
        M = throt(abs(rf_10(k)), angle(rf_10(k))) * M;
        M = A*M+B;
        M = zrot(phi)*M;
    end 
    Msig_fa10(j) = M(1) + 1i*M(2);
    Mss_fa10(:, j) = M;
end

figure
hold on; grid on;
plot(position*1e3, Mss_fa10(1, :))
plot(position*1e3, Mss_fa10(2, :))
plot(position*1e3, Mss_fa10(3, :))
plot(position*1e3, abs(Msig_fa10(:)))
legend('Mx', 'My', 'Mz', 'Mxy')
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization vs position FA=10 df=0')

% FT of the RF
rf_10_FT = fftshift(fft(rf_10(1:2001))); % the later half is unwanted
rf_30_FT = fftshift(fft(rf_30(1:2001)));
rf_90_FT = fftshift(fft(rf_90(1:2001)));

% compare the FT and slice profile
figure
hold on; grid on;
plot(position*1e3, abs(Msig_fa10(:)))
plot(position*1e3, abs(Msig_fa30(:)))
plot(position*1e3, abs(Msig_0(:))) % fa = 90
plot(position*1e3, abs(rf_10_FT(900:1100)))
plot(position*1e3, abs(rf_30_FT(900:1100)))
plot(position*1e3, abs(rf_90_FT(900:1100)))
legend('bloch 10', 'bloch 30', 'bloch 90', 'FT 10', 'FT 30', 'FT 90')
xlabel('position (mm)')
ylabel('Magnetization')
title("BlochSim - FT comparison T2=100ms")
% title("BlochSim - FT comparison T2=2ms")


%% Q3.4 

% Just simulate the time right after RF pulse ends ~ 2ms FA = 90 deg
fa = 1/2* pi;
rf_amplitude = fa/(2*pi*gamma*rf_area)*100;  % Scaling factor for the sinc
rf = rf_amplitude * rf_normalized;
Gz = Gss * ones(1, length(t));
df = 0;

Msig_noRephase = zeros(1, length(position));
Mss_noRephase = zeros(3, length(position)); % 0Hz df

% Simulation
for j = 1: length(position)
    M = [0, 0, 1]';
    [A, B] = freeprecess(1000*dt/2, T1, T2, df);
    for k = 1:length(rf)
        M = A*M + B;
        phi = 2*pi*gamma*position(j)*Gz(k)*dt/2;
        M = zrot(phi)*M;
        M = throt(abs(rf(k)), angle(rf(k))) * M;
        M = A*M+B;
        M = zrot(phi)*M;
    end 
    Msig_noRephase(j) = M(1) + 1i*M(2);
    Mss_noRephase(:, j) = M;
end

figure
hold on; grid on;
plot(position*1e3, Mss_noRephase(1, :))
plot(position*1e3, Mss_noRephase(2, :))
plot(position*1e3, Mss_noRephase(3, :))
plot(position*1e3, abs(Msig_noRephase))
legend('Mx', 'My', 'Mz', 'Mxy')
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization vs position FA=90 df=0 No rephasing')


%% Q3.5 SMS acquisition
% 25mm is the periodicity, BW corrsponding to 5mm
slice_offsets = [-2, -1, 0, 1, 2] * BW* 5; % frequency offset
rf_sms = zeros(1,length(rf));

for indx=1:length(slice_offsets)
    rf_sms = rf_sms + rf .* exp(2i*pi*slice_offsets(indx).*t);
end

rf_sms = [rf_sms 0*rf_sms];
Gz = [Gz -0.5*Gz];

Msig_SMS = zeros(1, length(position));
Mss_SMS = zeros(3, length(position)); % 0Hz df
position = [-6:0.01:6]*1e-2;% -60 to 60 mm
df = 0;

% Simulation
for j = 1: length(position)
    M = [0, 0, 1]';
    [A, B] = freeprecess(1000*dt/2, T1, T2, df);
    for k = 1:length(rf_sms)
        M = A*M + B;
        phi = 2*pi*gamma*position(j)*Gz(k)*2*dt/2;
        M = zrot(phi)*M;
        M = throt(abs(rf_sms(k)), angle(rf_sms(k))) * M;
        M = A*M+B;
        M = zrot(phi)*M;
    end 
    Msig_SMS(j) = M(1) + 1i*M(2);
    Mss_SMS(:, j) = M;
end

figure
hold on; grid on;
plot(position*1e3, Mss_SMS(1, :))
plot(position*1e3, Mss_SMS(2, :))
plot(position*1e3, Mss_SMS(3, :))
plot(position*1e3, abs(Msig_SMS))
legend('Mx', 'My', 'Mz', 'Mxy')
xlabel('position (mm)')
ylabel('Magnetization')
title('Magnetization vs position FA=90 df=0 SMS')