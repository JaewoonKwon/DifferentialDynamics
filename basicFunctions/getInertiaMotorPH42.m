function phi_PH42 = getInertiaMotorPH42()
% PH42-20-S300-R 
m_PH42 = 0.34;
com_PH42 = 1e-3 * [3.5664868e-02; -2.0225854e+00; -4.0415121e+01];
iner_PH42 = 1e-9 * [8.0053974e+05;7.9243138e+05;8.4106655e+04;-1.3999674e+02;6.5696853e+02;-3.8552358e+04];
phi = [m_PH42;zeros(3,1);iner_PH42];
SE3_com = eye(4);
SE3_com(1:3,4) = -com_PH42;
phi_PH42 = GToPhi_latest(largeAdjoint(SE3_com)' * PhiToG_latest(phi) * largeAdjoint(SE3_com));
end