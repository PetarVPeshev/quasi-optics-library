function [eta_tap, eff_area] = taper_efficiency(Efeed_comp, k, ...
    sph_grid, varargin)
%TAPER_EFFICIENCY Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(varargin{1}, 'reflector')
        focal_distance = varargin{2};
        diameter = varargin{3};
        physical_area = pi * (diameter / 2) ^ 2;
    
        [e_aperture, cyl_grid, ~] = aperture_field(Efeed_comp, k, ...
            sph_grid, 'reflector', focal_distance, diameter);
    
        drho = cyl_grid(1, 2, 1) - cyl_grid(1, 1, 1);
        dphi = cyl_grid(2, 1, 2) - cyl_grid(1, 1, 2);
    
        eff_area_nom = sum( sum( e_aperture .* cyl_grid(:, :, 1) ) ) ...
            * drho * dphi;
        eff_area_nom = abs( eff_area_nom(:, :, 1) ^ 2 ...
            + eff_area_nom(:, :, 2) ^ 2 + eff_area_nom(:, :, 3) ^ 2 );
        eff_area_denom = abs( e_aperture(:, :, 1) .^ 2 ...
            + e_aperture(:, :, 2) .^ 2 + e_aperture(:, :, 3) .^ 2 );
        eff_area_denom = sum( sum( eff_area_denom .* ...
            cyl_grid(:, :, 1) ) ) * drho * dphi;
        eff_area = eff_area_nom / eff_area_denom;
    
        eta_tap = eff_area / physical_area;
    end
end