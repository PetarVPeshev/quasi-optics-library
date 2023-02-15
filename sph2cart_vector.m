function cart_vector = sph2cart_vector(sph_vector, varargin)
%SPH2CART_VECTOR This function converts vector in spherical
%coordinate system to a vector in cartesian coordinate system
%   The function takes M-by-N-by-3 matrix for the spherical vector, and 
%   M-by-N matrix for Theta and Phi. The spherical matrix's third dimension
%   represents the spherical unit vectors as radial, elevation angle, and
%   azimuth angle vectors respectively. The M-by-N matricies of each
%   spherical unit vector holds its magnitude for the corresponding vector
%   location at the elevation and azimuth. The function returns M-by-N-by-3
%   matrix for the cartesian vector with third dimensions representing the
%   x, y, and z unit vectors respectively.
    if length(varargin) == 4
        for idx = 1:2:4
            if strcmp(varargin{idx}, 'Theta')
                theta = varargin{idx + 1};
            else
                phi = varargin{idx + 1};
            end
        end
    elseif length(varargin) == 2
        theta = varargin{1};
        phi = varargin{2};
    elseif length(varargin) == 1
        theta = varargin{1}(:, :, 1);
        phi = varargin{1}(:, :, 2);
    else
        error('Invalid arguments');
    end
    
    % Cartesian vector computation
    cart_vector = NaN( [size(sph_vector)] );
    cart_vector(:, :, 1) = ...
        sph_vector(:, :, 1) .* sin(theta) .* cos(phi) + ...
        sph_vector(:, :, 2) .* cos(theta) .* cos(phi) - ...
        sph_vector(:, :, 3) .* sin(phi);
    cart_vector(:, :, 2) = ...
        sph_vector(:, :, 1) .* sin(theta) .* sin(phi) + ...
        sph_vector(:, :, 2) .* cos(theta) .* sin(phi) + ...
        sph_vector(:, :, 3) .* cos(phi);
    cart_vector(:, :, 3) = ...
        sph_vector(:, :, 1) .* cos(theta) - ...
        sph_vector(:, :, 2) .* sin(theta);
end

