function uv_grid = uv_repr(varargin)
%UV_REPR This function calculates UV coordinate representation
%   This function takes the elevation and azimuth angle and calculates the
%   UV representation of the provided spherical coordinate points. For a 3D
%   matrix input, the third dimension represents the elevation and azimuth
%   angles respectively. For other input formats, the coordinates are
%   either specified or in the order elevation and azimuth angle. The
%   function returns the UV representation in a M-by-N-by-2 matrix with
%   third dimension representing U and V respectively.
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
    
    % UV coordinate computation
    uv_grid = NaN( [size(theta), 2] );
    uv_grid(:, :, 1) = sin(theta) .* cos(phi);
    uv_grid(:, :, 2) = sin(theta) .* sin(phi);
end