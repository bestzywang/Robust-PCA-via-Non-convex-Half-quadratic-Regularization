function [Out_X, RMSE,peaksnr,U,V,RRMSE,N] = RPCA_HQF(M_noise,rak,maxiter,ip,lam_1)
% This matlab code implements 'Robust PCA via Non-convex Half-quadratic Regularization'
%M - original data without noise
% M_noise - m x n observed matrix corrupted by noise
% rak is the rank of the object matrix
% maxiter - maximum number of iterations
% ip - controls the sparse
% If you have any questions, please contact Dr. Wang Zhiyong (z.y.w.@my.cityu.edu.hk)
[m,n] = size(M_noise);
U = rand(m,rak);
Ip = 3;
in = 0;
RMSE = [];
RRMSE = [];
peaksnr = [];
lam_2 = lam_1;
% PF find the initialization
for iter = 1 : Ip
    PINV_U = pinv(U,1e-8);
    V= PINV_U * M_noise;
   
    PIMV_V = pinv(V,1e-8);
    U = M_noise * PIMV_V;
end
X = U*V;
T = M_noise - X;    
t_m_n = T(:);
scale = 10*1.4815*median(abs(t_m_n - median(t_m_n)));
% scale = 10*std(t_m_n);
sigma = ones(m,n)*diag(scale);
ONE_1=ones(m,n);
ONE_1(abs(T-median(t_m_n))-sigma<0)=0;
N = T.*ONE_1;
U_p = U;
V_p = V;
for iter = 1 : maxiter
    D = M_noise - N;
    U = (D*V'-lam_1*U_p)*inv(V*V'-lam_1);
    V = inv(U'*U-lam_2)*(U'*D-lam_2*V_p);
    U_p = U;
    V_p = V;
    X = U*V;

    T = M_noise - X;    
    t_m_n = T(:);
    scale =min([scale ip*1.4815*median(abs(t_m_n - median(t_m_n)))]);

    sigma = ones(m,n)*scale;
    ONE_1=ones(m,n);
    ONE_1(abs(T-median(t_m_n))-sigma<0)=0;
    N = T.*ONE_1;
    
    RMSE= [RMSE norm(M_noise-X,'fro')/sqrt(m*n)];
    if iter~=1
        step_MSE = RMSE(iter-1) - RMSE(iter);
        if step_MSE < 0.000001
            in = in  + 1;
        end
        if in > 1
            break;
        end
    end
end
    Out_X = X;
end