function ad = smallAdjoint(se3)
if(size(se3,1)~=6 || size(se3,2) ~=1)
    error('[ERROR] Wrong input for smallAdjoint()')
end
w_bracket   = skew(se3(1:3));
v_bracket   = skew(se3(4:6));
ad          = [w_bracket, zeros(3, 3); v_bracket, w_bracket];
end

