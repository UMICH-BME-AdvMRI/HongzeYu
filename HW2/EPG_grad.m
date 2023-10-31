function [FpFmZ] = EPG_grad(FpFmZ, flag)
%   Function for gradient in EPG propagation
%   FpFmZ - EPG states [Fn, F-n, Zn]
%   flag - indicate the sign of the gradient

FpFmZ = [FpFmZ [0;0;0]];	% Add higher dephased state.

if (flag == 1)
    FpFmZ(1,:) = circshift(FpFmZ(1,:),[0 1]);	% Shift Fp states.
    FpFmZ(2,:) = circshift(FpFmZ(2,:),[0 -1]);	% Shift Fm states.
    FpFmZ(2,end)=0;					            % Zero highest Fm state.
    FpFmZ(1,1) = conj(FpFmZ(2,1));	            % Fill in lowest Fp state.
end

if (flag == -1)
    FpFmZ(2,:) = circshift(FpFmZ(2,:),[0 1]);	% Shift Fm states.
    FpFmZ(1,:) = circshift(FpFmZ(1,:),[0 -1]);	% Shift Fp states.
    FpFmZ(1,end)=0;					            % Zero highest Fp state.
    FpFmZ(2,1) = conj(FpFmZ(1,1));			    % Fill in lowest Fm state.
end
