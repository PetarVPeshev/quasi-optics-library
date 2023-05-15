function ibf = bf_current(relat_permit, dipole_l, dipole_w, dx, dy, k0, ...
    sph_grid_p, varargin)
%BF_CURRENT Summary of this function goes here
%   Detailed explanation goes here
    te_amplitude = 1;
    tm_amplitude = 1;
    if ~isempty(varargin)
        te_amplitude = varargin{1};
        tm_amplitude = varargin{2};
    end
    
    [k_comp, k] = wave_vector(relat_permit, k0, sph_grid_p);

    B = 2 * k * ...
        (cos(- k_comp(:, :, 1) * dipole_l / 2) - cos(k * dipole_l / 2)) ...
        .* sinc(- k_comp(:, :, 2) * dipole_w / (2 * pi)) ...
        ./ ( (k ^ 2 - (- k_comp(:, :, 1)) .^ 2) * sin(k * dipole_l / 2) );
    voltage = (- sin(sph_grid_p(:, :, 2)) * te_amplitude ...
        + cos(sph_grid_p(:, :, 2)) * tm_amplitude) .* B / sqrt(dx * dy);

    [km_comp, ~, ~] = floquet_modes(k, relat_permit, dx, dy, sph_grid_p);
    sgf = floquet_sgf(relat_permit, k, km_comp);
    impedance = array_impedance(dipole_l, dipole_w, dx, dy, k, km_comp, ...
        sgf);

    ibf = voltage ./ impedance;
end

%     k = k0 * sqrt(relat_permit);
%     k_comp = NaN(1, 1, 2);
%     k_comp(:, :, 1) = k * sin(sph_grid_p(:, :, 1)) ...
%         .* cos(sph_grid_p(:, :, 2));
%     k_comp(:, :, 2) = k * sin(sph_grid_p(:, :, 1)) ...
%         .* sin(sph_grid_p(:, :, 2));

%     incid_grid = NaN(1, 1, 2);
%     incid_grid(1, 1, 1) = theta_p;
%     incid_grid(1, 1, 2) = phi_p;

