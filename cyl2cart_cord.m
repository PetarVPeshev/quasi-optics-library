function cart_cord = cyl2cart_cord(varargin)
%CYL2CART_CORD Wrap function of pol2cart, which returns the cartesian
%coordinates in a 3D matrix
%   Detailed explanation goes here
    if length(varargin) == 6
        for idx = 1:2:6
            if strcmp(varargin{idx}, 'Rho')
                rho = varargin{idx + 1};
            elseif strcmp(varargin{idx}, 'Phi')
                phi = varargin{idx + 1};
            else
                z = varargin{idx + 1};
            end
        end
    elseif length(varargin) == 3
        rho = varargin{1};
        phi = varargin{2};
        z = varargin{3};
    elseif length(varargin) == 1
        rho = varargin{1}(:, :, 1);
        phi = varargin{1}(:, :, 2);
        z = varargin{1}(:, :, 3);
    else
        error('Invalid arguments');
    end
    
    % Cartesian coordinates computation
    cart_cord = NaN( [size(rho), 3] );
    [cart_cord(:, :, 1), cart_cord(:, :, 2), ...
        cart_cord(:, :, 3)] = pol2cart(phi, rho, z);
end

