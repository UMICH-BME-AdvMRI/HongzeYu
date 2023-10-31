function [m, sig] = EPG_single_echo_SE(m, T1, T2, TE, TR)
%   Function for simulating one TR of single echo SE using EPG

m = EPG_RF(m, pi/2, 0); % 90x excitation pulse

% Assume crusher gradient for 1 unit cycle twist
m = EPG_relax(m, T1, T2, TE/2); % relaxation before crusher t = TE/2
m = EPG_grad(m, 1); % crusher gradient
m = EPG_RF(m, pi, pi/2); % 180 pulse y-axis
m = EPG_grad(m, 1); % crusher gradient
m = EPG_relax(m, T1, T2, TE/2); % relaxation after crusher

sig = abs(m(1,1)); % signal acquired for this pixel
m = EPG_relax(m, T1, T2, TR - TE); % relaxation till the end of TR
