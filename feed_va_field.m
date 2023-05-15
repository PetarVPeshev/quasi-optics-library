function [va_field, rad_power] = feed_va_field(radius, k, sph_grid, ...
    varargin)
%FEED_VA_FIELD Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;
    dx = 0;
    kx = 0;
    if ~isempty(varargin)
        dx = varargin{1};
        kx = varargin{2};
    end

%     airy_pattern = 2 * (pi * (radius ^ 2)) ...
%         .* besselj(1, k * sin(sph_grid(:, :, 1)) * radius) ...
%         ./ (k * sin(sph_grid(:, :, 1)) * radius);

    Jeq = ft_current(k, radius, sph_grid(:, :, 1), 'circular', 'y');
    Jeq = Jeq / 2;

    va_field = zeros( [size(sph_grid, 1, 2), 3] );
    va_field(:, :, 2) = Jeq(:, :, 2) .* cos(sph_grid(:, :, 1)) ...
        .* sin(sph_grid(:, :, 2));
    va_field(:, :, 3) = Jeq(:, :, 2) .* cos(sph_grid(:, :, 2));

    va_field = - 1j * k * wave_impedance * va_field .* exp(1j * kx * dx);
    
%     dtheta = sph_grid(1, 2, 1) - sph_grid(1, 1, 1);
%     dphi = sph_grid(2, 1, 2) - sph_grid(1, 1, 2);

%     va_total_field = total_field(va_field);
%     rad_power = sum( sum( (va_total_field .^ 2) ...
%         .* sin(sph_grid(:, :, 1)) * dtheta * dphi ) ) ...
%         / (2 * wave_impedance);

    [~, ~, rad_power] = directivity(1, va_field, sph_grid, 1);
end

