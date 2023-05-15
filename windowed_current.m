function jw = windowed_current(dipole_l, dipole_w, i0, array_factor, k, ...
    k_comp)
%WINDOWED_CURRENT Summary of this function goes here
%   Detailed explanation goes here
    I = 2 * k * ( cos(k_comp(:, :, 1) * dipole_l / 2) ...
        - cos(k * dipole_l / 2) ) ...
        ./ ( (k ^ 2 - k_comp(:, :, 1) .^ 2) * sin(k * dipole_l / 2) );
    J = sinc(k_comp(:, :, 2) * dipole_w / (2 * pi));
    
    jw = zeros( [size(k_comp, 1, 2), 3] );
    jw(:, :, 1) = i0 .* I .* J .* array_factor(:, :, 1) ...
        .* array_factor(:, :, 2);
end

