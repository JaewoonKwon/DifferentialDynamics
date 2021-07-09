function phi_PH54 = getInertiaMotorPH54()
% PH54
m_PH54 = 0.855;
com_PH54 = 1e-3 * [2.1450651e+00;-3.2654927e+00;-5.7444049e+01];
iner_PH54 = 1e-9 * [4.0487283e+06;4.0350489e+06;3.4893573e+05;1.6444686e+04;1.4186275e+05;-2.1783575e+05];
phi = [m_PH54;zeros(3,1);iner_PH54];
SE3_com = eye(4);
SE3_com(1:3,4) = -com_PH54;
phi_PH54 = GToPhi_latest(largeAdjoint(SE3_com)' * PhiToG_latest(phi) * largeAdjoint(SE3_com));
end