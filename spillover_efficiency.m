function eta_so = spillover_efficiency(Efeed, relat_permit, sph_grid, ...
    r, varargin)
%SPILLOVER_EFFICIENCY Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(varargin{1}, 'reflector')
        focal_distance = varargin{2};
        diameter = varargin{3};
    
        theta_max = 2 * acot( 4 * focal_distance / diameter );
        [~, valid_col] = find(sph_grid(:, :, 1) <= theta_max);
        valid_col = unique(valid_col);
        
        Efeed_on_r = Efeed(:, valid_col, :);
        sph_grid_on_r = NaN( [size(sph_grid, 1), length(valid_col), 2] );
        sph_grid_on_r(:, :, 1) = sph_grid(:, valid_col, 1);
        sph_grid_on_r(:, :, 2) = sph_grid(:, valid_col, 2);
        dtheta_on_r = sph_grid_on_r(1, 2, 1) - sph_grid_on_r(1, 1, 1);
        dphi_on_r = sph_grid_on_r(2, 1, 2) - sph_grid_on_r(1, 1, 2);
    
        dtheta = sph_grid(1, 2, 1) - sph_grid(1, 1, 1);
        dphi = sph_grid(2, 1, 2) - sph_grid(1, 1, 2);
    
        [~, rad_intensity_on_r, ~] = directivity(relat_permit, ...
            Efeed_on_r, sph_grid_on_r, r);
        [~, rad_intensity, ~] = directivity(relat_permit, Efeed, ...
            sph_grid, r);
    
        eta_so_nom = sum( sum( rad_intensity_on_r ...
            .* sin(sph_grid_on_r(:, :, 1)) ) ) * dtheta_on_r * dphi_on_r;
        eta_so_denom = sum( sum( rad_intensity ...
            .* sin(sph_grid(:, :, 1)) ) ) * dtheta * dphi;
        eta_so = eta_so_nom / eta_so_denom;
    end
end