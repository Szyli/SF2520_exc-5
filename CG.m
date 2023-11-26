function [iter_num,rel_res_size] = CG(A,b,x_k,eps,K)
    %formula: lecture 15; page: 9/28 & 11/28
    %diagonal of A; in lecture notes it's named 'M'; is sparse
    D = diag(diag(A));
    T = -(A-D); %lecture 15; page: 3/28; is sparse
    x = x_k + 1;
    p_k = x_k;
    %alpha_k = 0;
    beta_k = 0;
    r_k = b;    %so we don't overwrite it in the iteration basically
    %iteration
    rel_res_size = [];
    iter_num = 0;
    while norm(r_k) > eps
        p_k = r_k + beta_k*p_k;
        alpha_k = (p_k'*r_k)/(p_k'*A*p_k);
        if iter_num > 0
            beta_k = (r_k'*r_k)/(r_k_m1'*r_k_m1);
        end
        r_k_m1 = r_k;   %r_k_m1 = r_{k-1}; from the formula
        r_k = r_k - alpha_k*(A*p_k);
        x_k = x_k + alpha_k * p_k;
        
        iter_num = iter_num+1;
        rel_res_size = [rel_res_size, normest(A*x_k-b)/norm(b)];
        if iter_num > K
            break
        end
    end
end

