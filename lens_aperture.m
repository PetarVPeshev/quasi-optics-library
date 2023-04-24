function [Ja, lens, cart_grid] = lens_aperture(diameter, f, n, ...
    lens_relat_permit, relat_permit, cyl_grid, varargin)
    %LENS_APERTURE Summary of this function goes here
    %   Detailed explanation goes here
    if isempty(varargin)
        lens = lens_parameters(diameter, lens_relat_permit, cyl_grid);
    elseif length(varargin) == 1
        lens = lens_parameters(diameter, lens_relat_permit, cyl_grid, ...
            varargin{1});
    else
        error('Error. Invalid number of arguments.');
    end

    wave_impedance = 376.730313668 / sqrt(relat_permit);

    S = sqrt(cos(lens.theta_t) .* (lens.eccint ...
        * cos(lens.sph_grid(:, :, 2)) - 1) ./ ( cos(lens.theta_i) ...
        .* (lens.eccint - cos(lens.sph_grid(:, :, 2))) ));
    
    [Ei, ~] = lens_feed(f, lens_relat_permit, lens.sph_grid(:, :, 1), ...
        lens.sph_grid(:, :, 2 : 3), n);

    [par_coeff, per_coeff] = transm_coeff(lens.theta_i, lens.theta_t, ...
        relat_permit, lens_relat_permit);

    Ja_const = - 2 * S ./ (wave_impedance * lens.sph_grid(:, :, 1));
    Ja = zeros( [size(Ja_const, 1, 2), 3] );
    Ja(:, :, 1) = Ja_const .* par_coeff .* Ei(:, :, 2);
    Ja(:, :, 2) = Ja_const .* per_coeff .* Ei(:, :, 3);
    Ja = cyl2cart_vector(Ja, lens.cyl_grid);

    cyl_grid(:, :, 3) = zeros( [size(cyl_grid, 1, 2)] );
    cart_grid = cyl2cart_cord(cyl_grid);
end

