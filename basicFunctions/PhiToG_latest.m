function G = PhiToG_latest(Phi)
% // Exception
if(size(Phi,1) ~= 10)
    error(['ERROR: Wrong input size for PhiToG(): Phi.rows() = ' num2str(size(Phi,1))])
    G = nan(6,6);
    return
end

% Phi = mass, h, ixx, iyy, izz, ixy, ixz, iyz
G = zeros(6,6);
m_Eye = Phi(1) * eye(3);
h_bracket = ToMatrix(Phi(2:4));
I_moment = [Phi(5), Phi(8), Phi(9);
            Phi(8), Phi(6), Phi(10);
            Phi(9), Phi(10), Phi(7)];

%     // substitute
G(1:3,1:3) = I_moment;
G(1:3,4:6) = h_bracket;
G(4:6,1:3) = h_bracket';
G(4:6,4:6) = m_Eye;
end