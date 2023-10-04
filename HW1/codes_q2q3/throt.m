function Rth=throt(phi,theta)
%   Function for rotation of angle phi around axis in transverse plane
%   defined by y=x*tan(theta)

Rth = zrot(theta)*xrot(phi)*zrot(-theta);