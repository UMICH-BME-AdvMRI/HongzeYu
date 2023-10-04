function Mss = Mss_bSSFP(fa,T1,T2,TE,TR,df)
%	Calculate the steady state magnetization at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  df is the resonant
%	frequency in Hz.  fa is flip angle in radians.

%   The codes are referencing http://mrsrl.stanford.edu/~brian/bloch/

Rflip = yrot(fa);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);
[Ate,Bte] = freeprecess(TE,T1,T2,df);

% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)

Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte); % SS magnetization