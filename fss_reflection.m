function [G, T] = fss_reflection(relat_permit, dipole_l, dipole_w, ...
    dx, dy, ibf, k, k_comp, sph_grid, varargin)
%FSS_REFLECTION Summary of this function goes here
%   Detailed explanation goes here
    te_amplitude = 1;
    tm_amplitude = 1;
    x = 0;
    y = 0;
    z = 0;
    for idx = 1 : 1 : length(varargin)
        if strcmp(varargin{idx}, 'Vte')
            te_amplitude = varargin{idx + 1};
        elseif strcmp(varargin{idx}, 'Vtm')
            tm_amplitude = varargin{idx + 1};
        elseif strcmp(varargin{idx}, 'x')
            x = varargin{idx + 1};
        elseif strcmp(varargin{idx}, 'y')
            y = varargin{idx + 1};
        elseif strcmp(varargin{idx}, 'z')
            z = varargin{idx + 1};
        end
    end
    
%     sgf = dyadic_sgf(relat_permit, k, k_comp, 'E', 'J');
    B = 2 * k * ...
        (cos(- k_comp(:, :, 1) * dipole_l / 2) - cos(k * dipole_l / 2)) ...
        .* sinc(- k_comp(:, :, 2) * dipole_w / (2 * pi)) ...
        ./ ( (k ^ 2 - (- k_comp(:, :, 1)) .^ 2) * sin(k * dipole_l / 2) );

    [km_comp, ~, ~] = floquet_modes(k, relat_permit, dx, dy, sph_grid);
    sgf = floquet_sgf(relat_permit, k, km_comp);
    
    e_inc = (- sin(sph_grid(:, :, 2)) * te_amplitude ...
        + cos(sph_grid(:, :, 2)) * tm_amplitude) .* exp(-1j ...
        * (k_comp(:, :, 1) .* x + k_comp(:, :, 2) .* y ...
        - k_comp(:, :, 3) .* z)) / sqrt(dx * dy);
    e_scat = ibf .* sgf(:, :, 1, 1) .* B .* exp(-1j ...
        * (km_comp(:, :, 1) .* x + km_comp(:, :, 2) .* y ...
        + km_comp(:, :, 3) .* z)) / (dx * dy);
    
    G = e_scat ./ e_inc;
    T = 1 + G;
end

