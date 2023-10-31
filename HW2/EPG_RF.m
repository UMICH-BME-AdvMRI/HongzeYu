function [FpFmZ] = EPG_RF(FpFmZ,alpha,phi)
%   Function for RF rotation in EPG propagation
%   FpFmZ - EPG states [Fn, F-n, Zn]
%   R - RF rotation matrix
%   alpha - Flip Angle in radians
%   phi - rotation axis from Mx in radians

R = [(cos(alpha/2))^2 exp(2i*phi)*(sin(alpha/2))^2 -1i*exp(1i*phi)*sin(alpha);
      exp(-2i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 1i*exp(-1i*phi)*sin(alpha);
      -1i/2*exp(-1i*phi)*sin(alpha) 1i/2*exp(1i*phi)*sin(alpha)      cos(alpha)];
FpFmZ = R * FpFmZ;
