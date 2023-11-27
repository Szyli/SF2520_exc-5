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
pause;
clear all;
%% Task 2

% Load the file by double-taping it, it will be loaded as A; otherwise
% it's a struct if you use load("...").

% spy(A);

% Don't forget to load the file
%% Task 2a)
NxN = size(A);
b = rand(NxN(1), 1);
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
semilogy(0:length(resVec)-1,resVec/norm(b))
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
[~, ~, ~, ~, resVec1] = pcg(A, b, 1e-4, 900, M1);
toc

L = ichol(A);
disp('Time 2')
tic
[x, flag, relRes, Iter, resVec2] = pcg(A, b, 1e-4, 900, L, L');
toc

semilogy(0:length(resVec1)-1,resVec1/norm(b))
hold on;
semilogy(0:length(resVec2)-1,resVec2/norm(b))
hold off;
xlabel('Iterations k')
grid;
ylabel('Relative residual')
legend('M = diag(A)', 'M = LL^T')
title('Convergence Scheme of `pcg`')

%% Solving the preconditioning equations
M = diag(diag(A));
M = inv(M);

disp('M = diag(A)')
tic
x1 = (M*A)\(M*b);
toc
disp('M = L*L^T')
L = ichol(A);
M = L*L';
M = inv(M);
tic
x2 = (M*A)\(M*b);
toc

%{
    M = diag(A);
        1.687100s (varies)
    M = L*L^T;
        0.794070s (varies)
%}

%% Task 2c)

%{
    Load the file convdiff.mat the same way as prior
%}
clear all;
%% pcg doesn't work
NxN = size(A);
b = rand(NxN(1), 1);

[x, flag, relRes, Iter, resVec] = pcg(A, b, 1e-4, 2000);
disp(flag)
% flag = varies inbetween 4 & 1; diverges for 4.

%% GMRES 
disp('Time for Gmres')
tic
[x, flag, relRes, Iter, resVec] = gmres(A, b, [], 1e-4, 500);
toc
semilogy(0:length(resVec)-1,resVec/norm(b))
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
[x, flag, relRes, Iter, resVec1] = gmres(A, b, [], 1e-4, 500, M);
toc
% Time: 5.191848s
%% 
disp('M = L*U factorization')
[L, U] = ilu(A);

tic
[x, flag, relRes, Iter, resVec2] = gmres(A, b, [], 1e-4, 500, L, U);
toc
disp(relRes)
% Time: 0.720025s

%% Residual error gmres with preconditioning.
semilogy(0:length(resVec1)-1,resVec1/norm(b))
hold on;
semilogy(0:length(resVec2)-1,resVec2/norm(b))
hold off;
xlabel('Iterations k')
grid;
ylabel('Relative residual')
legend('M = diag(A)', 'M = LU')
title('Convergence Scheme of `gmres`')

%% Using backslash

M = diag(diag(A));
M = inv(M);
tic
x1 = (M*A)\(M*b);
toc

%%
[L, U] = ilu(A);
M = L*U;
M = inv(M);
tic
x2 = (M*A)\(M*b);
toc
%{
    M = diag(A);
    : 3.635695s


    M = L*U;
    :  

%}

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
    

