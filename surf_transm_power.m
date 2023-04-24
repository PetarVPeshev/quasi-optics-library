function [te_power, tm_power] = surf_transm_power(par_coeff, per_coeff, ...
    theta_i, theta_t, med1_relat_permit, med2_relat_permit)
%SURF_TRANSM_POWER Summary of this function goes here
%   Detailed explanation goes here
    med_1_imped = 376.730313668 / sqrt(med1_relat_permit);
    med_2_imped = 376.730313668 / sqrt(med2_relat_permit);

    te_power = ( abs(per_coeff) .^ 2 ) .* med_2_imped .* cos(theta_t) ...
        ./ ( med_1_imped .* cos(theta_i) );
    tm_power = ( abs(par_coeff) .^ 2 ) .* med_2_imped .* cos(theta_t) ...
        ./ ( med_1_imped .* cos(theta_i) );
end

