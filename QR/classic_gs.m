function [Q, R] = clgs(A, print)
    [m, n]=size(A);
    Q = zeros(m,n);
    R = zeros(n,n);
     % the classical Gram-Schmidt
    for j = 1:n
        v = A(:, j);
        for i = 1:(j-1)
            R(i, j) = dot(Q(:, i), A(:, j));
            v = v - R(i, j)*Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
    if print
        disp(sprintf('||A - Q*R||_2 = %0.5e', norm(A - Q*R))); 
        disp(sprintf('||Q''*Q - I||_2 = %0.5e',norm(Q'*Q -A ))); % test orthogonality of Q
    end
end
