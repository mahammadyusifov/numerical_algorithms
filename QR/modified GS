function [Q, R] = mgs(A, print)
    [m, n]=size(A);
    Q = zeros(m,n);
    R = zeros(n,n);
    % the modified Gram-Schmidt
    V = zeros(size(A));
    for i = 1 : n
        V(:, i) = A(:, i);
    end
    for i = 1 : n
        R(i, i) = norm(V(:, i));
        Q(:, i) = V(:, i) / norm(R(i, i));
        for j = i + 1 : n
            R(i, j) = dot(Q(:, i), V(:, j));
            V(:, j) = V(:, j) - R(i, j) * Q(:, i);
        end
    end
    if print
        disp(sprintf('||A - Q*R||_2 = %0.5e', norm(A - Q*R))); 
        disp(sprintf('||Q''*Q - I||_2 = %0.5e', norm(Q'*Q - eye(n)))); % test orthogonality of Q
    end
end
