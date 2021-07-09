function adstar = smallAd_star(wrench)
if(length(wrench)~=6)
    error('ERROR: ad_star() must take as input a 6-dim wrench.');
    adstar = nan;
end

tau = wrench(1:3);
f   = wrench(4:6);
adstar = zeros(6,6);
adstar(1:3,1:3) = skew(tau);
adstar(1:3,4:6) = skew(f);
adstar(4:6,1:3) = skew(f);
end