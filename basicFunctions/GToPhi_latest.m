function Phi = GToPhi_latest(G)
% // Exception
if(size(G,1) ~= 6 || size(G,2) ~= 6)
    error('ERROR: Wrong input size for GToPhi_latest().')
end
% Phi = mass, h, ixx, iyy, izz, ixy, ixz, iyz
m = G(4,4);
h = ToVector(G(1:3,4:6));
I = [G(1,1); G(2,2); G(3,3); G(1,2); G(1,3); G(2,3)];

%     // substitute
Phi = [m; h; I];
end