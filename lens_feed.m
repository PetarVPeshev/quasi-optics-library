function [E, rad_power] = lens_feed(n, r, relat_permit, sph_grid, varargin)
%LENS_FEED Summary of this function goes here
%   Detailed explanation goes here
% FIXME: fix argument parsing
    const = ( cos(sph_grid(:, :, 1)) .^ n ) ./ r;
    if ~isempty(varargin)
        if ~strcmp(varargin{1}, 'NeglectPhase')
            k = sqrt(relat_permit) * 2 * pi * varargin{1} ...
                / physconst('LightSpeed');
            const = const .* exp(-1j * k .* r);
        end
    end

    E = zeros( [size(sph_grid, 1, 2), 3] );
    E(:, :, 2) = const .* cos(sph_grid(:, :, 2));
    E(:, :, 3) = - const .* sin(sph_grid(:, :, 2));
    
    [~, ~, rad_power] = directivity(relat_permit, E, sph_grid, r);
end
