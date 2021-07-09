function matrix = ToMatrix(vector)
if (length(vector) == 6)
    matrix = [0, -vector(3), vector(2), vector(4);
                vector(3), 0, -vector(1), vector(5);
                -vector(2), vector(1), 0, vector(6);
                0, 0, 0, 0];
    return
elseif (length(vector) == 3)
    matrix = [ 0, -vector(3), vector(2);
    vector(3), 0, -vector(1);
    -vector(2), vector(1), 0];
    return
else
    matrix = nan(4, 4);
    error(['ToMatrix(const VectorXd &vector) error: vector.rows() = ' num2str(size(vector,1))])
    return
end
end