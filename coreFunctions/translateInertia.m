function dPhi_vel = translateInertia(phi, spaceVel)

% % INPUTS
% phi = mass, CoM, inertia (in the following order: xx, yy, zz, xy, xz, yz) (10 x 1)
% spaceVel: se3 velocity of rigid body motion expressed in the world frame

% % OUTPUTS
% dPhi_vel: derivative of phi w.r.t. the rigid body motion given by spaceVel

if any(size(phi)~=[10 1]) || any(size(spaceVel)~=[6 1]) 
    error(['Wrong input: size(phi)=' num2str(size(phi)) ', size(spaceVel)=' num2str(size(spaceVel))])
end
m = phi(1);
h = phi(2:4);
i11 = phi(5);
i22 = phi(6);
i33 = phi(7);
i12 = phi(8);
i13 = phi(9);
i23 = phi(10);
w = spaceVel(1:3);
v = spaceVel(4:6);
dPhi_vel = zeros(10,1);
dPhi_vel(2:4) = -skew(h)*w + m*v;
hv = h.*v;
i12w3 = i12*w(3);
i13w2 = i13*w(2);
i23w1 = i23*w(1);
dPhi_vel(5) = -2*(i12w3 - i13w2 - hv(3) - hv(2)); 
dPhi_vel(6) = -2*(-i12w3 + i23w1 - hv(3) - hv(1)); 
dPhi_vel(7) = -2*(i13w2 - i23w1 - hv(2) - hv(1)); 
dPhi_vel(8) = i11*w(3) - i13*w(1) -i22*w(3) + i23*w(2) - h(2)*v(1) - h(1)*v(2); 
dPhi_vel(9) = -i11*w(2) + i12*w(1) -i23*w(3) + i33*w(2) - h(3)*v(1) - h(1)*v(3); 
dPhi_vel(10) = -i12*w(2) + i22*w(1) + i13*w(3) - i33*w(1) - h(3)*v(2) - h(2)*v(3); 
end