function [msig,m]=slice_profile(rf,grad,t,T1,T2,pos,df)

gamma = 4258;
dt = t(2)-t(1);         % s.
rfrot = 2*pi*gamma*rf*dT;
pos = pos(:).';		
msig = 0*pos;
m = [msig;msig;msig];

for x=1:length(pos)
    M = [0;0;1];
    [A,B] = freeprecess(1000*dT/2,T1,T2,df);

    for k = 1:length(rf)
	    M = A*M+B;
	    grot = 2*pi*gamma*(pos(x)/10)*grad(k)*dT/2;
	    M = zrot(grot)*M;
        M = throt(abs(rfrot(k)),angle(rfrot(k))) * M;
	    M = A*M+B;

	    grot = 2*pi*gamma*(pos(x)/10)*grad(k)*dT/2;
	    M = zrot(grot)*M;
    end
    m(:,x) = M;
    msig(x) = M(1)+1i*M(2);
end