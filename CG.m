function [x_k, iter_num, rel_res_size] = CG(A, b, x_k, eps, K)
    %formula: lecture 15; page: 9/28 & 11/28   

    p_k = b;
    
    beta_k = 0;

    r_k = b;    %so we don't overwrite it in the iteration basically
    
    rel_res_size = [];

    iter_num = 0;

    while norm(r_k) > eps && iter_num < K
        %%{
        A_pk = A * p_k;

        alpha_k = (p_k'*r_k)/(p_k'*A_pk);

        x_k = x_k + alpha_k * p_k;

        r_k = r_k - alpha_k*(A_pk);

        if iter_num > 0
            beta_k = (r_k'*r_k)/(r_k_m1'*r_k_m1);
        end

        r_k_m1 = r_k;

        p_k = r_k + beta_k*p_k;
        
        iter_num = iter_num + 1;

        %}
        
        %{
        p_k = r_k + beta_k*p_k;

        A_pk = A * p_k;

        alpha_k = (p_k'*r_k)/(p_k'*A_pk);

        if iter_num > 0
            beta_k = (r_k'*r_k)/(r_k_m1'*r_k_m1);
        end

        r_k_m1 = r_k;   %r_k_m1 = r_{k-1}; from the formula
        r_k = r_k - alpha_k*(A_pk);
        x_k = x_k + alpha_k * p_k;

        %}

        rel_res_size = [rel_res_size, norm(r_k)];
    end
end

