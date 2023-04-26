function [J, lens, cart_grid, varargout] = lens_aperture(diameter, n, ... 
    lens_relat_permit, relat_permit, cyl_grid, varargin)
    %LENS_APERTURE Summary of this function goes here
    %   Detailed explanation goes here
    lens = lens_parameters(diameter, lens_relat_permit, cyl_grid);
    enable_ft = false;
    if ~isempty(varargin)
        for idx = 1 : 1 : length(varargin)
            if strcmp(varargin{idx}, 'FT')
                enable_ft = true;
                k_comp = varargin{idx + 1};
            elseif isequal(size(varargin{idx}), [1 1])
                lens = lens_parameters(diameter, lens_relat_permit, ...
                    cyl_grid, varargin{idx});
            end
        end
    elseif length(varargin) > 3
        error('Error. Invalid number of arguments.');
    end

    wave_impedance = 376.730313668 / sqrt(relat_permit);

    S = sqrt(cos(lens.theta_t) .* (lens.eccint ...
        * cos(lens.sph_grid(:, :, 2)) - 1) ./ ( cos(lens.theta_i) ...
        .* (lens.eccint - cos(lens.sph_grid(:, :, 2))) ));
    
    [Ei, ~] = lens_feed(n, lens.sph_grid(:, :, 1), lens_relat_permit, ...
        lens.sph_grid(:, :, 2 : 3), 'NeglectPhase');

    [par_coeff, per_coeff] = transm_coeff(lens.theta_i, lens.theta_t, ...
        relat_permit, lens_relat_permit);

    Ja_const = - 2 * S ./ (wave_impedance);
    J = zeros( [size(Ja_const, 1, 2), 3] );
    J(:, :, 1) = Ja_const .* par_coeff .* Ei(:, :, 2);
    J(:, :, 2) = Ja_const .* per_coeff .* Ei(:, :, 3);
    J = cyl2cart_vector(J, lens.cyl_grid);

    cyl_grid(:, :, 3) = zeros( [size(cyl_grid, 1, 2)] );
    cart_grid = cyl2cart_cord(cyl_grid);

    if enable_ft == true
        rho = cyl_grid(:, :, 1);
        drho = rho(1, 2) - rho(1, 1);
        phi = cyl_grid(:, :, 2);
        dphi = phi(2, 1) - phi(1, 1);

        num_pages = floor( memory().MaxPossibleArrayBytes * 0.8 / (96 * size(J, 1) * size(J, 2)) );
        while rem(size(k_comp, 1) * size(k_comp, 2), num_pages) ~= 0
            num_pages = num_pages - 1;
        end
        num_k_elements = size(k_comp, 1) * size(k_comp, 2);
        num_iterations = num_k_elements / num_pages;
        
        Jft = zeros( [size(k_comp, 1, 2), 3] );
        for coord_idx = 1 : 1 : 2
            Jcoord_idx = NaN(num_k_elements, 1);
            for iteration_idx = 1 : 1 : num_iterations
                start_idx = (iteration_idx - 1) * num_pages + 1;
                end_idx = iteration_idx * num_pages;

                kx = k_comp(start_idx : end_idx);
                ky = k_comp(num_k_elements + start_idx : num_k_elements + end_idx);
                kx = permute(repmat(kx, [1 1 size(J, 1) size(J, 2)]), [3 4 2 1]);
                ky = permute(repmat(ky, [1 1 size(J, 1) size(J, 2)]), [3 4 2 1]);
                
                Jiteration = sum( sum( J(:, :, coord_idx) .* exp( 1j * kx .* rho .* cos(phi) ) .* exp( 1j * ky .* rho .* sin(phi) ) .* rho ) ) * drho * dphi;
                Jcoord_idx(start_idx : end_idx) = permute(Jiteration(1, 1, :), [3 1 2]);
            end
            Jft(:, :, coord_idx) = reshape(Jcoord_idx, [size(k_comp, 1) size(k_comp, 2)]);
        end

        varargout{1} = Jft;
    end
end

