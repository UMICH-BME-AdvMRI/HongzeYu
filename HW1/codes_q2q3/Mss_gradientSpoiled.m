function Mss = Mss_gradientSpoiled(fa,T1,T2,TE,TR,df,dephasing_moment,perfect_spoiler)
%	Calculate the steady state magnetization at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  df is the resonant
%	frequency in Hz.  fa is flip angle in radians.
%   Dephasing_moment represents the dephasing angle given by the spoiler

%   The codes are referencing http://mrsrl.stanford.edu/~brian/bloch/

Rflip = yrot(fa);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);
[Ate,Bte] = freeprecess(TE,T1,T2,df);

if perfect_spoiler
    perfect_spoiler = [0 0 0; 0 0 0; 0 0 1];
    Mss = inv(eye(3)-Ate*Rflip*Atr*perfect_spoiler) * (Ate*Rflip*Btr+Bte); 
else
    %   At the end of each TR, add dephasing given by spoiler
    Atr = zrot(dephasing_moment) * Atr;    
    Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte); % SS magnetization
end