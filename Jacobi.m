function [x_k, iter_num, rel_res_size] = Jacobi(A, b, x_k, eps, K)
    %diagonal of A; in lecture notes it's named 'M'; is sparse
    M = diag(diag(A));
    T = M - A;
    x_guess = x_k + 1;
    A_x = A * x_guess;
    

    rel_res_size = norm(A*x_k - b, 2)/norm(b,2);
    iter_num = 1;

    while norm(A_x - b,2)/norm(b) > eps && iter_num < K  %there is no do while in matlab...
        x_guess = x_k;
        %lecture 14; page: 18/20 
        x_k = M\(T*x_guess + b);
        iter_num = iter_num + 1;
        A_x = A * x_k;
        rel_res_size = [rel_res_size; norm(A_x - b, 2)/norm(b,2)];
    end
    
    %spectral_radius = max(eig(M\T)); disp(spectral_radius)
end
