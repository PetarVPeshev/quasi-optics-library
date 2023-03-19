function [a, b, c, eccint, max_angle] = lens_aperture(diameter, ...
    relat_permit)
%LENS_APERTURE Summary of this function goes here
%   Detailed explanation goes here
    eccint = 1 / sqrt(relat_permit);
    
    crit_angle = asin(eccint);
    max_angle = pi / 2 - crit_angle;
    
    r_min = diameter / (2 * sin(max_angle));
    
    a = r_min * (1 - eccint * cos(max_angle)) / (1 - eccint ^ 2);
    c = a * eccint;
    b = sqrt(a ^ 2 - c ^ 2);
end

