function V_reg = V_regressor_latest(V)
if(length(V)~=6)
    error('ERROR: V_regressor(V) takes only a 6-dim. vector as input!');
    V_reg = nan;
end
V_reg = zeros(6,10);
w = V(1:3);
v = V(4:6);

V_reg(4:6,1) = v;
V_reg(1:3,2:4) = -skew(v);
V_reg(4:6,2:4) = skew(w);
V_reg(1:3,5:7) = diag(w);
V_reg(1:3,8:10)= [
    w(2),w(3), 0;
    w(1), 0, w(3);
    0, w(1), w(2)
];

end