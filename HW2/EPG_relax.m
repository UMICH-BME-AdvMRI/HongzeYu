function FpFmZ = EPG_relax(FpFmZ, T1, T2, T)
%   Relaxation in EPG propagation
%   FpFmZ - EPG states [Fn, F-n, Zn]
%   T1,T2 - relaxation in ms
%   T - time in ms

% Relaxation
E = diag([exp(-T/T2) exp(-T/T2) exp(-T/T1)]);
FpFmZ = E * FpFmZ;
FpFmZ(3,1) = FpFmZ(3,1) + (1 - exp(-T/T1));