function phi_XH540 = getInertiaMotorXH540()
% XH540 
m_XH540 = 0.165;
com_XH540 = 1e-3 * [-4.0527847e-1;-1.9957815e1;-2.1294880e1];
iner_XH540 = 1e-9 * [7.3026397e+04;3.9200049e+04;6.0676644e+04;-1.4986434e+02;-2.0521503e+02;-8.7305765e+03];
phi = [m_XH540;zeros(3,1);iner_XH540];
SE3_com = eye(4);
SE3_com(1:3,4) = -com_XH540;
phi_XH540 = GToPhi_latest(largeAdjoint(SE3_com)' * PhiToG_latest(phi) * largeAdjoint(SE3_com));
end