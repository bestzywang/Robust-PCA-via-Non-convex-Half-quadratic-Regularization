clear
clc
close all
t_1 = [];
F_1_S = [];
RMSE_1 = [];
RMSE_1_S = [];
t_2 = [];
F_2_S = [];
RMSE_2 = [];
RMSE_2_S = [];
t_3 = [];
F_3_S = [];
RMSE_3 = [];
RMSE_3_S = [];
t_4 = [];
F_4_S = [];
RMSE_4 = [];
RMSE_4_S = [];
t_5 = [];
F_5_S = [];
RMSE_5 = [];
RMSE_5_S = [];
t_6 = [];
F_6_S = [];
RMSE_6 = [];
RMSE_6_S = [];

%%



for k = 1:100
   

    m = 500 ; rho_s = 0.1;
    n = m ;
    rak = m/50;
    U = (1*randn(m,rak)); V = (1*randn(n,rak));
    M = U*V' ;
    
    SNR = 6;
    noise = Gaussian_noise(M,'GM',SNR);
    M_noise = M + noise;
    S = noise;
    


    maxiter = 100;
    ip = 3;
    lam_1 = 0.001;
    tic
    [X_1, MSE,peak_snr,U,V,~,S_1] = RPCA_HQF(M_noise,rak,maxiter,ip,lam_1);
    toc
    t_1 = [t_1 toc];
    RMSE_1 = [RMSE_1 norm(X_1-M,'fro')/sqrt(m*n)];

end
mean_R_t = [mean(RMSE_1) std(RMSE_1) mean(t_1)]

