clear;

%initial data
%{
n = [10,20];  %sequence
d = [1,2];  %sequence

%}

n = [20,40,60];  %inner grid points; sequence
d = [1,2,3];  %dimensions; sequence
K = 1200;   %iteration limit

main_1(n,d,K,Methods.CG)


% iteration number is the length of data!!


% main function here
function main_1(n_seq,d_seq, K, method)
    if isinteger(K)
        disp('K must be an integer!')
        return
    end
    if isempty(n_seq) || isempty(d_seq)
        disp('n or d is empty sequence!')
        return
    end
    legends = [];
    if method == Methods.Jacobi
        for d = d_seq
            for n = n_seq
                %system terms
                N = n^d;    %should be >=2000
                b = rand(N,1);
                A = lap(N,1);   %already sparse
                x_k = zeros(N,1); %x_0 = 0
                eps = 0.001;

                titles = 'Jacobi method';
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
                
                semilogy(1:1:iter_num, rel_res_size)
                hold on
                legends = [legends; ['d=',num2str(d),' n=',num2str(n),'; ']];
            end
        end
    elseif method == Methods.CG
        for d = d_seq
            for n = n_seq
                %system terms
                N = n^d;    %should be >=2000
                b = rand(N,1);
                A = lap(N,1);   %already sparse
                x_k = zeros(N,1); %x_0 = 0
                eps = 0.001;

                titles = 'Conjugate gradient method';
                %{
                lambda = eig(A);
                kappa = abs(max(lambda))/abs(min(lambda));
                kappa = normest(A)*normest(inv(A)); %condition number of A
                %https://blogs.mathworks.com/cleve/2017/07/17/what-is-the-condition-number-of-a-matrix/
                beta = 1-(1/sqrt(kappa));
                disp(beta)
                %}
                
                [iter_num, rel_res_size] = CG(A,b,x_k,eps,K);
                
                semilogy(1:1:iter_num, rel_res_size)
                hold on
                legends = [legends; ['d=',num2str(d),' n=',num2str(n),'; ']];
            end
        end
    else
        disp('No valid method input; Methods.CG & Methods.Jacobi')
        return
    end
    ylabel('relative residual size')
    xlabel('iteration number')
    grid on
    hold off
    title(titles)
    legend(legends)
    pause;
    close;

end    
    

