function [ego, sph_grid] = go_field(pw_amplitude, k, sph_grid, angle_i, ...
    focal_distance, diameter)
%GO_FIELD Summary of this function goes here
%   Detailed explanation goes here
    theta_max = 2 * acot( 4 * focal_distance / diameter );
    [~, valid_col] = find(sph_grid(:, :, 1) <= theta_max);
    valid_col = unique(valid_col);
    sph_grid = sph_grid(:, valid_col, :);

    dkri = zeros(1, 3);
    dkri(1, 1) = k * sin(angle_i(1)) * cos(angle_i(2));
    dkri(1, 2) = k * sin(angle_i(1)) * sin(angle_i(2));

    kr = zeros( [size(sph_grid, 1, 2), 3] );
    kr(:, :, 1) = k * sin(sph_grid(:, :, 1)) .* cos(sph_grid(:, :, 2));
    kr(:, :, 2) = k * sin(sph_grid(:, :, 1)) .* sin(sph_grid(:, :, 2));

    delta = ( 1 - cos(sph_grid(:, :, 1)) ) ...
        ./ ( 1 + cos(sph_grid(:, :, 1)) );
    rho = focal_distance * dkri / k;
    rho = permute(repmat(rho, [1, 1, size(sph_grid, 1, 2)]), [3 4 2 1]);
    
    ego_zero_amplitude = - 2 * pw_amplitude ...
        ./ ( 1 + cos(sph_grid(:, :, 1)) );
    ego_zero = zeros( [size(sph_grid, 1, 2), 3] );
    ego_zero(:, :, 2) = ego_zero_amplitude .* sin(sph_grid(:, :, 2));
    ego_zero(:, :, 3) = ego_zero_amplitude .* cos(sph_grid(:, :, 2));

    ego = ego_zero .* exp( -1j * dot(kr, rho, 3) ) ...
        .* exp(-1j * dot(kr, rho, 3) .* delta);
end

