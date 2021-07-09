function A = getInertiaTransformRegressor(T)
R = T(1:3,1:3);
p = T(1:3,4);
if any(size(T)~=[4 4]) || any(any(R'*R-eye(3)>1e-9)) || det(R)-1>1e-9 || any(T(4,1:3) ~= zeros(1,3))
    error('Error 1')
end
A = zeros(10);
% mass
A(1,1) = 1;
% h
A(2:4,1) = p;
A(2:4,2:4) = R;
% I
idx = [1 5 9 4 7 8];
mReg = reshape(triu(-skew(p)^2),[],1);
A(5:10,1) = mReg(idx);
A(5:10,2:4) = -[0       -2*p(2) -2*p(3)
                -2*p(1) 0       -2*p(3)
                -2*p(1) -2*p(2) 0
                p(2)    p(1)    0
                p(3)    0       p(1)
                0       p(3)    p(2)] * R;
r1 = R(:,1);
r2 = R(:,2);
r3 = R(:,3);
IReg_xx = reshape(r1 * r1',[],1);
IReg_yy = reshape(r2 * r2',[],1);
IReg_zz = reshape(r3 * r3',[],1);
IReg_xy = reshape(r1 * r2' + r2 * r1',[],1);
IReg_xz = reshape(r1 * r3' + r3 * r1',[],1);
IReg_yz = reshape(r2 * r3' + r3 * r2',[],1);
A(5:10,5) = IReg_xx(idx);
A(5:10,6) = IReg_yy(idx);
A(5:10,7) = IReg_zz(idx);
A(5:10,8) = IReg_xy(idx);
A(5:10,9) = IReg_xz(idx);
A(5:10,10) = IReg_yz(idx);
end