function vec = ToVector(mat)
% se(3) matrix to se(3) vector (or so(3))

if size(mat,1) == 4          % se(3)
    vec = zeros(6,1);
    vec(1) = -mat(2,3);
    vec(2) = mat(1,3);
    vec(3) = -mat(1,2);
    vec(4) = mat(1,4);
    vec(5) = mat(2,4);
    vec(6) = mat(3,4);
elseif size(mat,1) == 3      % so(3)
    vec = zeros(3,1);
    vec(1) = -mat(2,3);
    vec(2) = mat(1,3);
    vec(3) = -mat(1,2);
else error('Check size of input matrix (It should be 3x3 or 4x4).')
    
end
