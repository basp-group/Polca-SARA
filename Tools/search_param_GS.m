%% Golden search


% DEFINE ALL YOUR PARAMETERS FOR YOUR ALGORITHM HERE
% param = ....
% x0 =  ...



% DEFINE YOUR ALGORITHM AS A FUNCTION
% par is the parameter you want to optimize
algoGS =@(par) - your_algo_GS(x0,param, par);
% -------------------------------------------------
% golden_search will minimize the first value sent by algoGS 
% -------------------------------------------------
% - if you want to maximize the SNR : 
%   algoGS =@(par) - your_algo_GS(x0,param, par);
%   the first variable sent by your_algo_GS must be the SNR
%   you have to put 'minus' since you want to MAXIMIZE the SNR
%   and golden_search will minimize the value sent by algoGS
% -------------------------------------------------
% - if you want to minimize the MSE : 
%   algoGS =@(par) your_algo_GS(x0,param, par);
%   the first variable sent by your_algo_GS must be the MSE
% -------------------------------------------------

nu_inf = 1e-8 ; % Lowerbound for GS
nu_sup = 1 ; % Upper bound for GS
[nu_opt, save_nu, save_SNR] = golden_search(algoGS, nu_inf, nu_sup) ;
nu = nu_opt ; % your optimized parameter

figure
subplot 211
plot(save_nu, 'o')
xlabel('iterations of Golden search algo')
ylabel('value of the parameter')
subplot 212
plot(save_SNR, 'x')
xlabel('iterations of Golden search algo')
ylabel('value of the corresponding SNR')

