function sph_cord = cart2sph_cord(varargin)
%CART2SPH_CORD Wrap function of cart2sph, which returns the spherical
%coordinates in a 3D matrix
%   The function takes the coordinate points separately or in M-by-N-by-3
%   matrix. For a M-by-N-by-3 matrix input, the third dimension represents 
%   x, y, and z coordinates respectively. For other input formats, the 
%   coordinates are either specified or in the order x, y, and z. The 
%   function returns the spherical coordinates in a M-by-N-by-3 matrix with
%   third dimension representing radial distance, elevation, and azimuth 
%   angle respectively.
    if length(varargin) == 6
        for idx = 1:2:6
            if strcmp(varargin{idx}, 'X')
                x = varargin{idx + 1};
            elseif strcmp(varargin{idx}, 'Y')
                y = varargin{idx + 1};
            else
                z = varargin{idx + 1};
            end
        end
    elseif length(varargin) == 3
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
    elseif length(varargin) == 1
        x = varargin{1}(:, :, 1);
        y = varargin{1}(:, :, 2);
        z = varargin{1}(:, :, 3);
    else
        error('Invalid arguments');
    end
    
    % Spherical coordinates computation
    sph_cord = NaN( [size(x), 3] );
    [sph_cord(:, :, 1), sph_cord(:, :, 2), ...
        sph_cord(:, :, 3)] = cart2sph(x, y, z);
end