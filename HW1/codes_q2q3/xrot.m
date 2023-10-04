function Rx=xrot(phi)
%   Rotation around x-axis

Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];