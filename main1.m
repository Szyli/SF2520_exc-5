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
[x3, y3] = main_1(200,2,K,Methods.CG);
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
%[x2, y2] = main_1(4000,1,K,Methods.Jacobi);
%[x3, y3] = main_1(200,2,K,Methods.Jacobi);
%[x4, y4] = main_1(20,3,K,Methods.Jacobi);

disp(y1)

%% Plotting
ylabel('relative residual size')
xlabel('iteration number')
grid on
semilogy(1:x1, y1)
hold on
%semilogy(1:x2, y2)
%semilogy(1:x3, y3)
%semilogy(1:x4, y4)
title('CG method')
%legend('2000:1', '4000:1', '200:2', '20:3');
ylabel('Reside r')
xlabel('Iteration K')
hold off;




% main function here
function [iter_num, rel_res_size] = main_1(n, d, K, method)
    if isinteger(K)
        disp('K must be an integer!')
        return
    end
    legends = [];
    if method == Methods.Jacobi

        %system terms
        N = n^d;    %should be >=2000
        b = rand(N,1);
        A = sparse(lap(n,d));   %already sparse
        x_k = zeros(N,1); %x_0 = 0
        eps = 0.001;

        %{
        D = zeros(size(A));
        if issymmetric(A)
            for k = 1:length(A)
                D(k,k) = A(k,k);
            end
        end
        D = sparse(D);
        %2-norm; normest is for sparse matrices
        %I don't know what beta is for yet
        beta = normest(D\(A-D));  %norm(inv(D)*(A-D));
        %}
                
        [iter_num, rel_res_size] = Jacobi(A,b,x_k,eps,K);
        
    elseif method == Methods.CG

        %system terms
        N = n^d;    %should be >=2000
        b = rand(N,1);
        A = sparse(lap(n,d));   %already sparse
        x_k = zeros(N,1); %x_0 = 0
        eps = 0.001;

        %{
        lambda = eig(A);
        kappa = abs(max(lambda))/abs(min(lambda));
        kappa = normest(A)*normest(inv(A)); %condition number of A
        %https://blogs.mathworks.com/cleve/2017/07/17/what-is-the-condition-number-of-a-matrix/
        beta = 1-(1/sqrt(kappa));
        disp(beta)
        %}
        
        [iter_num, rel_res_size] = CG(A,b,x_k,eps,K);
        
    else
        disp('No valid method input; Methods.CG & Methods.Jacobi')
        return
    end

end    
    

