function [e_aperture, cyl_grid, theta_max] = aperture_field(Efeed_comp, ...
    k, sph_grid, varargin)
%APERTURE_FIELD Summary of this function goes here
%   Detailed explanation goes here
    for idx = 1 : 1 : length(varargin)
        if strcmp(varargin{idx}, 'reflector')
            antenna = 'reflector';
            focal_distance = varargin{idx + 1};
            diameter = varargin{idx + 2};
        end
    end
    
    if strcmp(antenna, 'reflector')
        theta_max = 2 * acot( 4 * focal_distance / diameter );
        [~, valid_col] = find(sph_grid(:, :, 1) <= theta_max);
        valid_col = unique(valid_col);
    
        theta = sph_grid(:, valid_col, 1);
        cyl_grid = zeros( [size(theta, 1, 2), 3] );
        cyl_grid(:, :, 1) = 2 * focal_distance * tan(theta / 2);
        cyl_grid(:, :, 2) = sph_grid(:, valid_col, 2);
    
        r_prime = focal_distance ./ ( cos(theta / 2) .^ 2 );
    
        e_aperture = zeros( [size(cyl_grid, 1, 2), 3] );
        e_aperture(:, :, 1) = - Efeed_comp(:, valid_col, 2);
        e_aperture(:, :, 2) = - Efeed_comp(:, valid_col, 3);
        e_aperture = e_aperture .* exp(-2j * k * focal_distance) ...
            ./ r_prime;
        e_aperture = cyl2cart_vector(e_aperture, cyl_grid(:, :, 2));
    end
end

