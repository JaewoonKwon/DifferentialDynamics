function Ad = largeAdjoint(MatSE3)
if( any(size(MatSE3)~=[4,4]) )
    error('[ERROR] Wrong input for AdjointMatrix()')
end
R = MatSE3(1:3,1:3);
p = MatSE3(1:3,4);
pR = skew(p)*R;
Ad = [R, zeros(3, 3); pR, R];
end

