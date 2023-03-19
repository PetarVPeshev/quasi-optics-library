function [E, rad_power] = lens_feed(f, relat_permit, r, sph_grid, n)
%LENS_FEED Summary of this function goes here
%   Detailed explanation goes here
    c = physconst('LightSpeed');

    k = sqrt(relat_permit) * 2 * pi * f / c;
    const = ( cos(sph_grid(:, :, 1)) .^ n ) * exp(-1j * k * r) / r;

    E = zeros( [size(sph_grid, 1, 2), 3] );
    E(:, :, 2) = const .* cos(sph_grid(:, :, 2));
    E(:, :, 3) = - const .* sin(sph_grid(:, :, 2));

    rad_power = total_field(E) .^ 2;
end

