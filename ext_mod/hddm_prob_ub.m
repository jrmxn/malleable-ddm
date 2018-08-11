function op = hddm_prob_ub(v,a,z)
%"""Probability of hitting upper boundary."""
if v==0
    %otherwise op = nan since division by 0
    op = 0.5;
else
    op = (exp(-2*a*z*v) - 1) / (exp(-2*a*v) - 1);
end
end