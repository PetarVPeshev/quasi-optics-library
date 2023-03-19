function cart_cord = sph2cart_cord(varargin)
%SPH2CART_CORD Wrap function of sph2cart, which returns the spherical
%coordinates in a 3D matrix
%   This function takes the coordinate points separately or in M-by-N-by-3
%   matrix. For a M-by-N-by-3 matrix input, the third dimension reperesents
%   the radial distance, elevation angle, and azimuth angle respectively. 
%   For other input formats, the coordinates are either specified or in the
%   order radial distance, elevation angle, and azimuth angle. The function
%   returns the cartesian coordinates in a M-by-N-by-3 matrix with the
%   third dimension representing x, y, and z respectively.
    if length(varargin) == 6
        for idx = 1:2:6
            if strcmp(varargin{idx}, 'R')
                r = varargin{idx + 1};
            elseif strcmp(varargin{idx}, 'Theta')
                theta = varargin{idx + 1};
            else
                phi = varargin{idx + 1};
            end
        end
    elseif length(varargin) == 3
        r = varargin{1};
        theta = varargin{2};
        phi = varargin{3};
    elseif length(varargin) == 1
        r = varargin{1}(:, :, 1);
        theta = varargin{1}(:, :, 2);
        phi = varargin{1}(:, :, 3);
    else
        error('Invalid arguments');
    end
    
    % Cartesian coordinates computation
    cart_cord = NaN( [size(r), 3] );
    [cart_cord(:, :, 1), cart_cord(:, :, 2), ...
        cart_cord(:, :, 3)] = sph2cart(phi, theta, r);
end