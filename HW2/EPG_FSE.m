function [m, sig] = EPG_FSE(m, T1, T2, TR, ESP, ETL)
%   Function for simulating one TR of FSE using EPG
sig = zeros(1,ETL);
m = EPG_RF(m, pi/2, 0); % 90x excitation pulse

for echo_num=1:ETL
    % Assume crusher gradient for 1 unit cycle twist
    m = EPG_relax(m, T1, T2, ESP/2); % relaxation before crusher t = ESP/2
    m = EPG_grad(m, 1); % crusher gradient
    m = EPG_RF(m, pi, pi/2); % 180 pulse y-axis
    m = EPG_grad(m, 1); % crusher gradient
    m = EPG_relax(m, T1, T2, ESP/2); % relaxation after crusher
    sig(echo_num) = abs(m(1,1));
end
m = EPG_relax(m, T1, T2, TR - ESP*ETL);