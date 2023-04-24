function [ J, varargout ] = aperture_current(Efeed_comp, k, ...
    relat_permit, sph_grid, varargin)
%APERTURE_CURRENT_FT Summary of this function goes here
%   Detailed explanation goes here
% TODO: implement lens aperture here
    wave_impedance = 376.730313668 / sqrt(relat_permit);
    enable_ft = false;
    for idx = 1 : 1 : length(varargin)
        if strcmp(varargin{idx}, 'FT')
            enable_ft = true;
            k_comp = varargin{idx + 1};
        elseif strcmp(varargin{idx}, 'reflector')
            antenna = 'reflector';
            focal_distance = varargin{idx + 1};
            diameter = varargin{idx + 2};
        elseif strcmp(varargin{idx}, 'Love')
            formulation = 'Love';
        elseif strcmp(varargin{idx}, 'Schelkunoff')
            formulation = 'Schelkunoff';
        end
    end

    if strcmp(antenna, 'reflector')
        theta_max = 2 * acot( 4 * focal_distance / diameter );
        [~, valid_col] = find(sph_grid(:, :, 1) <= theta_max);
        valid_col = unique(valid_col);

        theta = sph_grid(:, valid_col, 1);
        phi = sph_grid(:, valid_col, 2);
        dphi = phi(2, 1) - phi(1, 1);
        rho = 2 * focal_distance * tan(theta / 2);
        drho = rho(1, 2) - rho(1, 1);

        r_prime = focal_distance ./ ( cos(theta / 2) .^ 2 );

        n = zeros(size(rho));
        n(:, :, 3) = 1;

        e_aperture = zeros( [size(rho, 1, 2), 3] );
        e_aperture(:, :, 1) = - Efeed_comp(:, valid_col, 2);
        e_aperture(:, :, 2) = - Efeed_comp(:, valid_col, 3);
        e_aperture = e_aperture .* exp(-2j * k * focal_distance) ...
            ./ r_prime;
        e_aperture = cyl2cart_vector(e_aperture, phi);  % Check here

        M = cross(e_aperture, n, 3);
        J = cross(n, cross(n, e_aperture, 3), 3) / wave_impedance;

        if strcmp(formulation, 'Schelkunoff')
            M = zeros( [size(J)] );
            J = 2 * J;
        elseif ~strcmp(formulation, 'Love')
            error('Error. Invalid argument.');
        end

        % NOTE: Current density Fourier Transform calculation is
        % implemented only for the electric current density, due to this
        % reason, it is performed only when Schelkunoff's formulation is
        % used; when the calculation of the magnetic current density
        % Fourier Transform is implemented, it should be set to be
        % performed for Love's formulation as well.
        % TODO: Implement aperture magnetic current density Fourier
        % Transform calculation
        % NOTE: The algorithm must be optimized further. The permute
        % function takes a long time, and the use of this function should 
        % be reduced. Further, the allocated memory requirement should be
        % reduced to increase the number of calculations in one iteration.
        % TODO: Optimize further the current Fourier Transform algorithm
        % NOTE: For number of elements in one dimension of the current
        % density than available memory, the algorithm can not allocate
        % number of pages, resulting in no calculations
        % FIXME: Implement a fix for number of elements in one dimension of
        % the current density larger than available memory bug by 
        % separating the elements of the current density into sections
        if enable_ft == true && strcmp(formulation, 'Schelkunoff')
%             tic
%             clear('formulation', 'focal_distance', 'diameter', 'theta_max', 'valid_col', 'theta', 'r_prime', 'n', 'e_aperture', 'Efeed_comp', 'k', 'relat_permit', 'sph_grid', 'varargin', 'wave_impedance');
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
%                     clear('kx', 'ky', 'Jiteration');
                end
                Jft(:, :, coord_idx) = reshape(Jcoord_idx, [size(k_comp, 1) size(k_comp, 2)]);
            end
%             toc
        end

        cart_grid = cyl2cart_cord( rho, phi, zeros(size(rho)) );
%         J = cyl2cart_vector(J, phi);
%         M = cyl2cart_vector(M, phi);

        if strcmp(formulation, 'Love')
            varargout{1} = M;
            varargout{2} = cart_grid;
        else
            varargout{1} = cart_grid;
        end

        if enable_ft == true && strcmp(formulation, 'Schelkunoff')
            varargout{length(varargout) + 1} = Jft;
        end
    elseif strcmp(antenna, 'lens')
        error('Not implemented');
    else
        error('Error. Invalid arguments.');
    end
end

