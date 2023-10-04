function Rz=zrot(phi)
%   Rotation around z-axis

Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
