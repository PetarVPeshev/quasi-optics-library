function [ego, vgo, sph_grid] = go_field(focal_distance, diameter, ...
    pw_amplitude, k, sph_grid, sph_grid_i)
%GO_FIELD Summary of this function goes here
%   Detailed explanation goes here
    theta_max = 2 * acot( 4 * focal_distance / diameter );
    [~, valid_col] = find(sph_grid(:, :, 1) <= theta_max);
    valid_col = unique(valid_col);
    sph_grid = sph_grid(:, valid_col, :);

    theta = sph_grid(:, :, 1);
    phi = sph_grid(:, :, 2);

    theta_i = sph_grid_i(:, :, 1);
    phi_i = sph_grid_i(:, :, 2);

    delta = ( 1 - cos(theta) ) ./ ( 1 + cos(theta) );

    rho_fp_x = repmat(focal_distance * sin(theta_i) .* cos(phi_i), ...
        [1, 1, size(sph_grid, 1, 2)]);
    rho_fp_x = permute(rho_fp_x, [3 4 1 2]);
    rho_fp_y = repmat(focal_distance * sin(theta_i) .* sin(phi_i), ...
        [1, 1, size(sph_grid, 1, 2)]);
    rho_fp_y = permute(rho_fp_y, [3 4 1 2]);

    krho = NaN( [size(theta, 1, 2), 2] );
    krho(:, :, 1) = k * sin(theta) .* cos(phi);
    krho(:, :, 2) = k * sin(theta) .* sin(phi);

    ego_zero_const = repmat(- 2 * pw_amplitude ./ (1 + cos(theta)), ...
        [1, 1, size(sph_grid_i, 1, 2)]);
    ego_zero = zeros( [size(sph_grid, 1, 2), size(sph_grid_i, 1, 2), 3] );
    ego_zero(:, :, :, :, 2) = ego_zero_const .* sin(phi);
    ego_zero(:, :, :, :, 3) = ego_zero_const .* cos(phi);

    dot_prod = krho(:, :, 1) .* rho_fp_x + krho(:, :, 2) .* rho_fp_y;
    ego = ego_zero .* exp(-1j * dot_prod .* (1 + delta));
    ego = permute(ego, [1 2 5 3 4]);
    vgo = ego * focal_distance / exp(1j * k * focal_distance);
end

%     dkri = zeros(1, 3);
%     dkri(1, 1) = k * sin(angle_i(1)) * cos(angle_i(2));
%     dkri(1, 2) = k * sin(angle_i(1)) * sin(angle_i(2));
% 
%     kr = zeros( [size(sph_grid, 1, 2), 3] );
%     kr(:, :, 1) = k * sin(sph_grid(:, :, 1)) .* cos(sph_grid(:, :, 2));
%     kr(:, :, 2) = k * sin(sph_grid(:, :, 1)) .* sin(sph_grid(:, :, 2));
% 
%     delta = ( 1 - cos(sph_grid(:, :, 1)) ) ...
%         ./ ( 1 + cos(sph_grid(:, :, 1)) );
%     rho = focal_distance * dkri / k;
%     rho = permute(repmat(rho, [1, 1, size(sph_grid, 1, 2)]), [3 4 2 1]);
%     
%     ego_zero_amplitude = - 2 * pw_amplitude ...
%         ./ ( 1 + cos(sph_grid(:, :, 1)) );
%     ego_zero = zeros( [size(sph_grid, 1, 2), 3] );
%     ego_zero(:, :, 2) = ego_zero_amplitude .* sin(sph_grid(:, :, 2));
%     ego_zero(:, :, 3) = ego_zero_amplitude .* cos(sph_grid(:, :, 2));
% 
%     ego = ego_zero .* exp( -1j * dot(kr, rho, 3) ) ...
%         .* exp(-1j * dot(kr, rho, 3) .* delta);
% 
%     vgo = ego * focal_distance / exp(1j * k * focal_distance);

