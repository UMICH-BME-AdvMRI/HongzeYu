%% Q2
% The codes are referencing http://mrsrl.stanford.edu/~brian/bloch/

%% Q2a Simulate steady-state frequency response of a bSSFP sequence

T1 = 1000;	% ms. Adopting from Q2b
T2 = 100;	% ms
TE = [2.5 5 10];	% ms.
TR = 2*TE;	% ms.
fa = pi/3;

df = [-100:100]; 	%Hz

Sig = zeros(length(df),length(TE));

for n=1:length(TE)

    for k=1:length(df)
		Mss = Mss_bSSFP(fa,T1,T2,TE(n),TR(n),df(k));
		Sig(k,n)=Mss(1)+1i*Mss(2);
    end

end

% Plots
figure
subplot(2,1,1);
plot(df,abs(Sig));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('TE=2.5', 'TE=5.0', 'TE=10');
title('q2a Steady State Frequency Response of bSSFP sequence')

subplot(2,1,2);
plot(df,angle(Sig));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis([min(df) max(df) -pi pi]);
grid on;
legend('TE=2.5', 'TE=5.0', 'TE=10');



%% q2b FLASH sequence

TE = 5;    % ms.
TR = 10;    % ms.

% q2b(i)
fa = 10/180 * pi;   % 10 degree in radians

Mss_perfect = Mss_gradientSpoiled(fa,T1,T2,TE,TR,0,0,1);
Sig_perfect=Mss_perfect(1)+1i*Mss_perfect(2);

% Display the results for Msig
fprintf('---Q2b(i)---\n')
fprintf('Msig Magnitude: %f\n', abs(Sig_perfect));
fprintf('Msig Phase: %f in radians\n', angle(Sig_perfect));

% Display the components of Mss
fprintf('Mss (Mx component): %f\n', Mss_perfect(1));
fprintf('Mss (My component): %f\n', Mss_perfect(2));
fprintf('Mss (Mz component): %f\n', Mss_perfect(3));


% q2b(ii)
N = 200;    % Number of spins to simulate
dephasing_moment = [2*pi 4*pi 8*pi 16*pi];
M_avg = zeros(3, length(dephasing_moment)); 
Sig_avg = zeros(length(dephasing_moment), 1); 

for n = 1:length(dephasing_moment)
    phi = ([1:N]/N-0.5) * dephasing_moment(n);
    M = zeros(3,N);     % Magnetization across voxel
    for k = 1:N
        Mss_gradSpoil = Mss_gradientSpoiled(fa,T1,T2,TE,TR,0,phi(k),0);
        M(:,k) = Mss_gradSpoil;
    end
    M_avg(:,n) = mean(M, 2); 
    Sig_avg(n) = abs(M_avg(1,n) + 1i*M_avg(2,n));
end

% Scatter the signal magnitude as a function of dephasing moment
figure;
scatter(dephasing_moment, round(Sig_avg,3), 'o');
xlabel('Dephasing Moment (radians)');
ylabel('Average Signal Magnitude');
title('q2b(ii) Signal Magnitude vs. Dephasing Moment');
ylim([0.0 1.0]);
grid off;
xticks(dephasing_moment); % Set specific tick locations for clarity
xticklabels({'2\pi','4\pi','8\pi','16\pi'}); % Label the ticks in terms of pi


% q2b(iii) RF Spoiling
dephasing_moment = 2*pi; % From the last question, 2*pi is sufficient
RF_phase_range = 0:1:180; % degrees, from 0 to 180 with a step of 1 degree
RF_phase_rad = deg2rad(RF_phase_range);
M_transverse = zeros(length(RF_phase_range), 1);
Nex = 200; % 200 excitations
df = 0;

for idx = 1:length(RF_phase_range)
    % Using the spgrsignal function for the RF spoiling simulation:
    [Msig, ~] = spgrsignal(fa, T1, T2, TE, TR, df, Nex, RF_phase_rad(idx));
    M_transverse(idx) = abs(Msig); % Store magnitude of transverse magnetization
end

% Scatter the transverse magnetization as a function of RF phase:
figure;
scatter(RF_phase_range, M_transverse, 'o');
xlabel('RF Phase (degrees)');
ylabel('Transverse Magnetization Magnitude');
title('q2b(iii)Transverse Magnetization vs. RF Phase');
grid on;

% Assuming M_transverse is already computed for each RF_phase_range_deg

[min_mag, min_idx] = min(M_transverse);
best_RF_phase_deg = RF_phase_range(min_idx);

fprintf('---Q2b(iii)---\n')
fprintf('The RF phase that best eliminates transverse magnetization is %f degrees.\n', best_RF_phase_deg);




