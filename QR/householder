function [Q, R] = house(A, print)
    [m, n]=size(A);
    V = zeros(m,n); % store reflection vectors
    % Householder 
    % generate Q using an implicit calculation of a product Qx
    Q = eye(m);
    B = A;
    for k = 1 : n
        e1 = zeros(m,1);
        e1(k) = 1;
        x = A(k : m, k);
        v = sign(x(1))*norm(x)*e1(k : m, 1) + x;
        v = v / norm(v);
        V(k : m, k) = v;
        A(k : m, k : n) = A(k : m, k : n) - 2 * v * v' * A(k : m, k : n);
        for t = 1 : m
            Q(:, t) = Q(:, t) - 2 * V(:,k) * (V(:, k)' * Q(:, t));
        end
    end
    R = A;
    Q = Q';
    if print
        disp(sprintf('||A - Q*R||_2 = %0.5e', norm(B - Q*R))); 
        disp(sprintf('||Q''*Q - I||_2 = %0.5e', norm(Q' * Q - eye(m)))); % test orthogonality of Q
    end
end
