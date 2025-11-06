function [angle_out] = pi22pi (angle_rad);

angle_rad = pi2pi (angle_rad);

f= find (angle_rad < 0);

angle_rad (f) = angle_rad(f) + (2*pi);
angle_out = angle_rad;
