function [te_power, tm_power] = surf_transm_power(par_coeff, per_coeff, ...
    theta_inc, relat_permit)
%SURF_TRANSM_POWER Summary of this function goes here
%   Detailed explanation goes here
    med_1_imped = 376.730313668;
    med_2_imped = 376.730313668 / sqrt(relat_permit);

    theta_tra = asin( sqrt(relat_permit) * sin(theta_inc) );

    te_power = ( abs(per_coeff) .^ 2 ) .* med_2_imped .* cos(theta_tra) ./ ...
        ( med_1_imped .* cos(theta_inc) );
    tm_power = ( abs(par_coeff) .^ 2 ) .* med_2_imped .* cos(theta_tra) ./ ...
        ( med_1_imped .* cos(theta_inc) );
end

