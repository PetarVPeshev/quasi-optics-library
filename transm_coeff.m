function [par_coeff, per_coeff, theta_tra] = transm_coeff(theta_inc, ...
    relat_permit)
%TRANSMISSION_COEFF Summary of this function goes here
%   Detailed explanation goes here
    med_1_imped = 376.730313668;
    med_2_imped = 376.730313668 / sqrt(relat_permit);

    theta_tra = asin( sqrt(relat_permit) * sin(theta_inc) );

    per_coeff = 2 * med_1_imped * cos(theta_inc) ./ ...
        (med_1_imped * cos(theta_inc) + med_2_imped * cos(theta_tra));
    par_coeff = 2 * med_1_imped * cos(theta_inc) ./ ...
        (med_1_imped * cos(theta_tra) + med_2_imped * cos(theta_inc));
end

