clear;

%initial data
%{
n = [10,20];  %sequence
d = [1,2];  %sequence

%}
%% CG
K = 1000;   %iteration limit
% [x, y] = main_1(n,d,K,Methods.CG)

[x1, y1] = main_1(2000,1,K,Methods.CG);
[x2, y2] = main_1(4000,1,K,Methods.CG);
[x3, y3] = main_1(50,2,K,Methods.CG);
[x4, y4] = main_1(20,3,K,Methods.CG);


%% Plotting
ylabel('relative residual size')
xlabel('iteration number')
grid on
semilogy(1:x1, y1)
hold on
semilogy(1:x2, y2)
semilogy(1:x3, y3)
semilogy(1:x4, y4)
title('CG method')
legend('2000:1', '4000:1', '200:2', '20:3');
ylabel('Reside r')
xlabel('Iteration K')
hold off;

% iteration number is the length of data!!
%% Jacobi

K = 1000;   %iteration limit
% [x, y] = main_1(n,d,K,Methods.CG)

[x1, y1] = main_1(2000,1,K,Methods.Jacobi);
[x2, y2] = main_1(4000,1,K,Methods.Jacobi);
[x3, y3] = main_1(200,2,K,Methods.Jacobi);
[x4, y4] = main_1(20,3,K,Methods.Jacobi);

%disp(y1);

%% Plotting
ylabel('relative residual size')
xlabel('iteration number')
semilogy(1:x1, y1)
hold on
semilogy(1:x2, y2)
semilogy(1:x3, y3)
semilogy(1:x4, y4)
title('Jacobi method')
grid on
legend('2000:1', '4000:1', '200:2', '20:3');
ylabel('Reside r')
xlabel('Iteration K')
hold off;

%% Task 2






%% Main function task 1.
function [iter_num, rel_res_size] = main_1(n, d, K, method)
    if isinteger(K)
        disp('K must be an integer!')
        return
    end

    N = n^d;    %should be >=2000
    b = rand(N,1);
    A = sparse(lap(n,d));   %already sparse
    x_k = zeros(N,1); %x_0 = 0
    eps = 1e-6;

    if method == Methods.Jacobi
                
        [iter_num, rel_res_size] = Jacobi(A, b, x_k , eps, K);
        
    elseif method == Methods.CG
        
        [iter_num, rel_res_size] = CG(A, b, x_k, eps, K);
        
    else
        disp('No valid method input; Methods.CG & Methods.Jacobi')
        return
    end

end    
    

