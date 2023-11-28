function [iter_num,rel_res_size] = Jacobi(A,b,x_k,eps,K)
    %diagonal of A; in lecture notes it's named 'M'; is sparse
    M = diag(diag(A));
    T = A - M;
    x_guess = x_k + 1;
    
    %iteration
    rel_res_size = norm(b-A*x_k,2)/norm(b,2);
    iter_num = 0;
    while norm(x_k - x_guess,2) > eps   %there is no do while in matlab...
        x_guess = x_k;
        %lecture 14; page: 18/20
        x_k = M\(b - T*x_k); %inv(D)*(T*x(:,k)+b); 

        iter_num = iter_num + 1;
        rel_res_size = [rel_res_size; norm(b-A*x_k,2)/norm(b,2)];
        
        if iter_num > K
            break
        end
    end
end

