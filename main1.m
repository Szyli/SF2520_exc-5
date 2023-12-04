clear;
%% CG
eps = 1e-6;
K = 2000;   %iteration limit
% [x, y] = main_1(n,d,K,Methods.CG)

[~, iter1, rel1] = main_1(40000,1,K,Methods.CG, eps);
[~, iter2, rel2] = main_1(200,2,K,Methods.CG, eps);
[~, iter3, rel3] = main_1(35,3,K,Methods.CG, eps);


%% Plotting
ylabel('relative residual size')
xlabel('iteration number')
grid on
semilogy(1:iter1, rel1)
hold on
semilogy(1:iter2, rel2)
semilogy(1:iter3, rel3)
title('CG method')
legend('40000:1', '200:2', '35:3');
ylabel('Reside r')
xlabel('Iteration K')
hold off;

%% Jacobi
eps = 1e-6;
K = 2000;   %iteration limit
% [x, y] = main_1(n,d,K,Methods.CG)

[~, iter1, rel1] = main_1(40000,1,K,Methods.Jacobi, eps);
[~, iter2, rel2] = main_1(200,2,K,Methods.Jacobi, eps);
[~, iter3, rel3] = main_1(35,3,K,Methods.Jacobi, eps);

%% Plotting
ylabel('relative residual size')
xlabel('iteration number')
semilogy(1:iter1, rel1)
hold on
semilogy(1:iter2, rel2)
semilogy(1:iter3, rel3)
title('Jacobi method')
grid on
legend('40000:1', '200:2', '35:3');
ylabel('Reside r')
xlabel('Iteration K')
hold off;
%% Questions
%{
    Both methods converges faster with smaller n values, and is more or
    less independent on the dimension d, in the aspect of numbers of
    iterations.

    Conjugate gradient method converges faster than the Jacobi method.

    When looking at the Jacobi method, we can identify that beta, beta =
    1-sigma, has some proportinally to the spectral radii, and thus when
    the spectral radii is greater or equal than one, then the solution does not
    converge. With smaller values of 1-sigma, fewer iterations are needed.
    For the conjugate gradient method, we define (1-1/sqrt(kappa)), where
    kappa is the condition number of the matrix, i.e. kappa = condest(A).
    This factor increases with larger matrices A, and thus more iterations
    are needed.

    With increasing d, with fixed N, the convergence is faster for both
    methods, however cg is faster.

%}

%%
clear all;

%% part b
clear;
eps = 1e-10;
K = 10000;   %iteration limit
n = 50;
d = 2;
disp('Cg')
tic
[x, iter, rel] = main_1(n,d,K,Methods.CG,eps);
toc

A = lap(n,d);
b = rand(n^d,1);
disp('\')
tic
x_exact = A\b;
toc
semilogy(1:iter, rel)
grid on
title('CG method')
legend('35:3');
ylabel('Reside r')
xlabel('Iteration K')

diff = norm(x - x_exact);
disp(diff)

%{
    CG is faster when n = 35, d = 3, time = 0.226s
    whilst \ time = 0.548s

    CGi is slower when n = 50, d = 2, time = 0.0142s
    whilst \ time = 0.00278s

    In order to check how the \ solves this system we use:
        spparms('spumoni', 2)
        A;
    From here it's visible that the solver uses a Cholesky solver. 
    From this, we theorize that the Cholesky solver solves the system:
        x = (L*L^(*))^(-1) * b.
    
    This then results in a time-complexity O(N * p^2) where p is the
    bandwith.

    If A is however banded, with some band-density a banded solver is
    instead used, which again has the time-complexity O(N * p^2).

    If we instead look at iterative methods, focusing on the conjugate
    gradient method; if the condition number is good enough, the number of
    iterations will be around C * (1-1/sqrt(kappa))^k, where k are the
    number of iterations. C is independent on k and kappa, thus in the case
    of kappa being good enough, the number of iterations are less than N *
    p^2, which then exceeds the 'stationary' method.   

%}
%% Task 2
clear;
% Load the file by double-taping it, it will be loaded as A; otherwise
% it's a struct if you use load("...").

% spy(A);

load cooling_flange.mat
NxN = size(A);
b = rand(NxN(1), 1);
spy(A);
%% Task 2a)

tic
[x, flag, relRes, Iter, resVec] = pcg(A, b, 1e-4, 900);
toc
disp('Converged')
disp(flag)
disp('Relative residual size')
disp(relRes)
disp('Iterations')
disp(Iter)
disp('Residual error size')
disp(size(resVec))
semilogy(0:Iter, resVec/norm(b))
xlabel('Iterations k')
grid;
ylabel('Relative residual')
title('Convergence Scheme of `pcg`')

% Computation time 0.878271.
%{
    Iterations;
        763 (Varies depending on b)
    Residual size;
        9.9645*10^(-5) (Varies depending on b)

%}
%% Solving with the inverse method
tic
x = A\b;
toc

% Computation time 1.679235.
% The computation time is more or less doubled.
%% Precondtioner diag(A) 2b)
M1 = diag(diag(A));

disp('Time 1')
tic
[~, ~, ~, iter1, resVec1] = pcg(A, b, 1e-4, 900, M1);
toc

L = ichol(A);
disp('Time 2')
tic
[x, flag, relRes, iter2, resVec2] = pcg(A, b, 1e-4, 900, L, L');
toc

semilogy(0:iter1,resVec1/norm(b))
hold on;
semilogy(0:iter2,resVec2/norm(b))
hold off;
xlabel('Iterations k')
grid;
ylabel('Relative residual')
legend('M = diag(A)', 'M = LL^T')
title('Convergence Scheme of `pcg`')

%% Solving the preconditioning equations
M = diag(diag(A));

disp('M = diag(A)')
tic
x1 = (M\A)\(M\b);
toc
%%
L = ichol(A);
disp("L*L' precond.")
tic
%x2 = (L'\L)\b;
toc

%{
    M = diag(A);
        1.514123s (varies)
    M = L*L^T;
        1732.55228s (varies) (29 min)
    
    However if we do L = inv(L) this is drastically faster.
%} 

%% Task 2c)

%{
    Load the file convdiff.mat the same way as prior
%}
clear all;
load convdiff.mat
NxN = size(A);
b = rand(NxN(1), 1);
%% pcg doesn't work
[~, flag, relRes, Iter, resVec] = pcg(A, b, 1e-4, 2000);
disp(flag)
% flag = varies inbetween 4 & 1; diverges for 4.

%% GMRES 
disp('Time for Gmres')
tic
[x, flag, relRes, Iter, resVec] = gmres(A, b, [], 1e-4, 500);
toc
semilogy(0:Iter(2), resVec/norm(b))
xlabel('Iterations k')
grid;
ylabel('Relative residual')
title('Convergence Scheme of `gmres`')

disp(relRes)
disp(Iter)
disp(size(resVec))

%{
    Time for gmres: 6.449367
    Relative residual size 9.9719 * 10^(-5)
    Iterations: 274
%}


%% Timing
disp('Time for \')
tic
x = A\b;
toc
% Time ~ 3.596066s

%% Using Preconditioners with gmres

M = diag(diag(A));
disp('M = diag(A)')
tic
[x, flag, relRes, Iter1, resVec1] = gmres(A, b, [], 1e-4, 500, M);
toc
% Time: 4.781603s
%% 
disp('M = L*U factorization')
[L, U] = ilu(A);

tic
[x, flag, relRes, Iter2, resVec2] = gmres(A, b, [], 1e-4, 500, L, U);
toc
disp(relRes)
% Time: 0.589242s

%% Residual error gmres with preconditioning.
semilogy(0:Iter1(2), resVec1/norm(M\b))
hold on;
semilogy(0:Iter2(2), resVec2/norm((L*U)\b))
hold off;
xlabel('Iterations k')
grid;
ylabel('Relative residual')
legend('M = diag(A)', 'M = LU')
title('Convergence Scheme of `gmres`')

%% Using backslash

M = diag(diag(A));
disp('Diag A \')
tic
x1 = (M\A)\(M\b);
toc

%%
[L, U] = ilu(A);
disp('LU \')
tic
% x2 = U\L\b;
toc
%{
    M = diag(A);
    : 3.635695s


    M = L*U;
    : undef (Ran out of ram)

%}

%% Main function task 1.
function [x, iter_num, rel_res_size] = main_1(n, d, K, method, eps)

    N = n^d;    %should be >=2000
    b = rand(N,1);
    A = sparse(lap(n,d));   %already sparse
    x_k = zeros(N,1); %x_0 = 0

    %disp('K-iterations')

    if method == Methods.Jacobi
                
        [x, iter_num, rel_res_size] = Jacobi(A, b, x_k , eps, K);
        
    elseif method == Methods.CG
        
        [x, iter_num, rel_res_size] = CG(A, b, x_k, eps, K);
        
    else
        disp('No valid method input; Methods.CG & Methods.Jacobi')
        return
    end

end    
    

