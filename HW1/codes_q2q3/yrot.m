function Ry=yrot(phi)
%   Rotation around y-axis

Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];