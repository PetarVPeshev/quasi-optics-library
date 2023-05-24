function lens = lens_parameters(diameter, relat_permit, cyl_grid, varargin)
%LENS_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    lens.diameter = diameter;
    lens.eccint = 1 / sqrt(relat_permit);

    if isempty(varargin)
        crit_angle = asin(lens.eccint);
        lens.theta_max = pi / 2 - crit_angle - 0.01 * pi / 180;
    elseif length(varargin) == 1
        lens.theta_max = varargin{1};
    else
        error('Error. Invalid number of arguments.');
    end
    
    lens.r_min = diameter / (2 * sin(lens.theta_max));
    
    lens.a = lens.r_min * (1 - lens.eccint * cos(lens.theta_max)) ...
        / (1 - lens.eccint ^ 2);
    lens.c = lens.a * lens.eccint;
    lens.b = sqrt(lens.a ^ 2 - lens.c ^ 2);

    lens.cyl_grid = NaN( [size(cyl_grid, 1, 2), 3] );
    lens.cyl_grid(:, :, 1 : 2) = cyl_grid;
    lens.cyl_grid(:, :, 3) = lens.c ...
        + lens.a * sqrt( 1 - (cyl_grid(:, :, 1) / lens.b) .^ 2 );

    lens.sph_grid = NaN( [size(lens.cyl_grid, 1, 2, 3)] );
    lens.sph_grid(:, :, 2) = atan(lens.cyl_grid(:, :, 1) ...
        ./ lens.cyl_grid(:, :, 3));
    lens.sph_grid(:, :, 1) = lens.a * (1 - lens.eccint ^ 2) ...
        ./ ( 1 - lens.eccint * cos(lens.sph_grid(:, :, 2)) );
    lens.sph_grid(:, :, 3) = lens.cyl_grid(:, :, 2);
    
    lens.theta_i = acos( ( 1 - lens.eccint ...
        .* cos(lens.sph_grid(:, :, 2)) ) ./ sqrt( 1 + lens.eccint .^ 2 ...
        - 2 * lens.eccint * cos(lens.sph_grid(:, :, 2)) ) );
    lens.theta_t = asin(sqrt(relat_permit) * sin(lens.theta_i));
end

