function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)
%
%	Function calculates the signal from an RF-spoiled sequence
%	after Nex excitations.
%
%   The codes are referencing http://mrsrl.stanford.edu/~brian/bloch/

Nf = 200;	% Simulate 200 different gradient-spoiled spins.
phi = [1:Nf]/Nf*2*pi;

M=zeros(3,Nf,Nex+1);

	
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);

M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);
	
Rfph = 0;
Rfinc = inc;

for n=1:Nex
    
	A = Ate * throt(flip,Rfph); % RF spoiling
	B = Bte;
	M = A*M+B*on;

	Msig = mean(squeeze(M(1,:)+1i*M(2,:))) * exp(-1i*Rfph); % demodulate the RF phase
	Mss = M;

	M=Atr*M+Btr*on;
    
    % Gradient spoiling
    for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
    end

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end
