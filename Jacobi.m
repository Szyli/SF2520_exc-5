function [iter_num,rel_res_size] = Jacobi(A,b,x_k,eps,K)
    %diagonal of A; in lecture notes it's named 'M'; is sparse
    M = diag(diag(A));
    
    T = A - M;
    
    x = x_k + 1;
    
    %iteration
    rel_res_size = [];
    iter_num = 0;
    while norm(b - A*x_k) > eps * norm(b)   %there is no do while in matlab...

        %lecture 14; page: 18/20
        x_k = M\(T*x_k+b); %inv(D)*(T*x(:,k)+b); 
        iter_num = iter_num + 1;
        %rel_res_size = [rel_res_size; norm(b-A*x_k)/norm(b)];
        rel_res_size = [rel_res_size; norm(b - x_k)];
        if iter_num > K
            break
        end
    end
end

