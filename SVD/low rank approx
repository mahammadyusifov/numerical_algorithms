function [Ak, k, rel_error] = lr_approx(A, epsilon, print)
    [U,S,V] = svd(A);
    k = 0;
    Sk= zeros(size(S));
    n = rank(A);
    sigma1 = norm(A);
    for i = 1 : n
        if S(i, i) > sigma1*epsilon
            Sk(i, i) = S(i, i);
            k = k + 1;
        end
    end
    Ak = U*Sk*V';
    rel_error = norm(S - Sk) / norm(S);
    if print
        disp(sprintf('k = %d', k)); 
        disp(sprintf('||A - A_k||_2/||A||_2 = %0.5e', rel_error)); 
    end
end
