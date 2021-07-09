function invSE3 = inverseSE3(SE3)
if(any(size(SE3)~=[4 4]))
    error('ERROR: Wrong input for InverseSE3()')
end
invSE3          = eye(4);
R               = SE3(1:3,1:3)';
invSE3(1:3,1:3) = R;
invSE3(1:3,4)   = - R * SE3(1:3,4);
end
