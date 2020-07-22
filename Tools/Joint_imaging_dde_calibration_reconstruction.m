function results = Joint_imaging_dde_calibration_reconstruction...
    (param_im, param_dde, param_algo, param_data)


% *************************************************************************
%
% Reconstruction algorithm for joint DDE calibration and Stokes imaging in RI
% Algorithm based on a block-coordinate forward-algorithm
%
% Related paper:
% A non-convex optimization algorithm for joint DDE calibration and imaging
% in radio interferometry
% Audrey Repetti, Jasleen Birdy, Arwa Dabbech, Yves Wiaux
% Arxiv preprint: 1701.03689
%
% Contact: a.repetti@hw.ac.uk
%
% *************************************************************************
%
% Same algorithm to estimate
% DIEs + image      if param_algo.initialization = 1
% or DDEs + image   if param_algo.initialization = 0
%
% *************************************************************************
% Output: structure -results- containing
%         eps_rec: reconstructed unknown sources
%         im_rec: reconstructed global image (x = xo + eps)
%         U1_rec/U2_rec: reconstructed DIEs or DDEs
%         
% *************************************************************************



%% Initialization


for i = 1:param_im.P
if param_algo.initialization == 1 % reconstruction with DIEs
xo{i} = param_im.xo{i} ;
epsilon{i} = 0*xo{i} ;
else      % reconstruction with DDEs
xo{i} = param_im.xo{i} ;
epsilon{i} = param_im.im_rec_die{i} - param_im.xo{i} ;
end
param_im.x{i} = xo{i} + epsilon{i} ;
Ntot = prod(param_im.Nt) ;
end

if strcmp(param_algo.opt_dict, 'Dirac')
prox_func =@(x,Ax) min(max(proxl1(x, param_algo.eta*(1.9/Ax) ),param_im.min_x{i}),Inf) ;
regul_x =@(x) param_algo.eta * sum( abs(x(:)) ) ;
else
param_l1.min_x = param_im.min_x ;
param_l1.max_x = param_im.max_x ;
param_l1.mask_app = (xo>0) ;
param_l1.Psi = param_algo.Psi ;
param_l1.Psit = param_algo.Psit ;
param_l1.real = 1 ;
param_l1.pos = 1 ;
param_l1.verbose = 0 ;
param_l1.weights = 1 ;
param_l1.xA = xo ;
prox_func =@(x,Ax) solver_prox_L1_full_image(x, param_algo.eta*1.9/Ax, param_l1) ;
regul_x =@(x) param_algo.eta * sum( abs(param_l1.Psit(xo+x)) ) ;
end
% ======================================================================= %




% -------------------------------------------------------------------------
% Some useful variables...
% -------------------------------------------------------------------------
param_dde.U1 = param_dde.U1_init ;
param_dde.U2 = param_dde.U2_init ;

D1 = cell(param_data.T,1) ;
D2 = cell(param_data.T,1) ;
YY1 = cell(param_data.T,1) ;
YY2 = cell(param_data.T,1) ;

for i =1:param_im.P
for s = 1:param_data.T
D1{s} = cell2mat(param_dde.U1{s}) ;
D2{s} = cell2mat(param_dde.U2{s}) ;
YY_temp1 = cell(1,param_data.na) ;
YY_temp2 = cell(1,param_data.na) ;
Y_noisy_temp{i} = param_data.Y{s,i} ;
parfor alpha = 1:param_data.na
YY_temp1{alpha} = transpose(Y_noisy_temp{i}(alpha, :)) ;
YY_temp2{alpha} = Y_noisy_temp{i}(:,alpha) ;
end
YY1{s,i} = YY_temp1 ; % observation to minimize wrt U1
YY2{s,i} = YY_temp2 ; % observation to minimize wrt U2
end
end
% -------------------------------------------------------------------------




SNR =@(x, xtrue) 20 * log10( sqrt( sum( xtrue(:).^2 ) / sum( (x(:)-xtrue(:)).^2 ) ) ) ;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Displaying variables
Time_tot = 0 ;
G = create_matrix_G_D1_D2_toolbox(param_data, param_dde, param_im)  ;
B_tmp   = @(x)   G * param_im.TF(x);  
Bt_tmp  = @(x)   real(param_im.TF_adj(G' * x));


for i =1:param_im.P
data_fid = 0.5* sum( abs( B_tmp(param_im.x{i}) - param_data.y{i}).^2 ) ;
reg_x = regul_x(epsilon{i}) ; 
reg_d = 0 ;
for s = 1:param_data.T
dk1_temp = param_dde.U1{s} ;
dk2_temp = param_dde.U2{s} ;
err_temp = cell(1,param_data.na) ;
parfor alpha = 1:param_data.na
err_temp{alpha} = sum( (abs(dk1_temp{alpha} - dk2_temp{alpha})).^2 ) ;
end
reg_d = reg_d + sum( cell2mat(err_temp) ) ;
end
crit = data_fid + param_dde.nu * reg_d ;
crit = crit + reg_x ;
error_dirty{i} = sqrt(sum( ( reshape(Bt_tmp(B_tmp(param_im.x{i}) - param_data.y {i}),prod(param_im.Ni-1),1)).^2 )) ;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% *************************************************************************
% *************************************************************************
% DISPLAY
% *************************************************************************
% *************************************************************************
disp(' ')
disp('Initialisation')
disp(['Crit: ',num2str(crit(1))])
disp(['Total time: ',num2str(sum(Time_tot))])
% disp(['Error dirty image: ',num2str((error_dirty(1)))])
disp(['SNR dirty image: ',num2str(SNR(Bt_tmp(B_tmp(param_im.x{i})), Bt_tmp(param_data.y{i})))])
disp(['SNR rec image: ',num2str(SNR(param_im.x{i}-param_im.xo{i}, param_im.im_true{i} - param_im.xo{i}))])
disp(' ')
% *************************************************************************
% *************************************************************************

end


%% Algorithm



for nit_tot = 1:param_dde.max_it
     
tic



%% ************************************************************************
% *************************************************************************
% Update DDEs    
% *************************************************************************
% *************************************************************************

% ------------------------------------------------------------------------
% create convolution matrix X
% -------------------------------------------------------------------------

for i = 1:param_im.P
X{i} = generate_convol_matrix_XX_toolbox(param_im, param_dde,i) ;
end
% -------------------------------------------------------------------------
    
for nbit_d = 1:param_dde.JUtot
U1_old = param_dde.U1 ;
U2_old = param_dde.U2 ;   

% -------------------------------------------------------------------------
% Update U1
% -------------------------------------------------------------------------

for i = 1:param_im.P

for s = 1:param_data.T
d2_temp = cell(1,param_data.na) ;
Supp_temp2 = param_dde.Support_Dalpha2{s} ;
dk2_temp = param_dde.U2{s} ;
for alpha = 1:param_data.na
d2_temp{alpha} = sparse(Ntot,1) ;
d2_temp{alpha}(Supp_temp2{alpha}) = dk2_temp{alpha}' ;
end
H_temp = X{i}*cell2mat(d2_temp) ;
YY_temp = YY1{s,i} ;
Supp_temp1 = param_dde.Support_Dalpha1{s} ;
dk1_temp = param_dde.U1{s} ;
parfor alpha = 1:param_data.na
% ----------------------------------------------------
% construction of the observation matrix 
H1_temp = transpose(H_temp(Supp_temp1{alpha}, :)) ;
H1_temp(alpha,:) = [] ;
% construction of data vector
YY_temp_alpha = YY_temp{alpha} ;
YY_temp_alpha(alpha,:) = [] ;
% initialization of u1
u1_temp = transpose(dk1_temp{alpha}(param_dde.Ind_conj)) ; 
% construction of sigma(u2) for regularization term
u2_temp = transpose(dk2_temp{alpha}(param_dde.Ind_conj)) ; 
% step size
Lips_temp = param_dde.nu + pow_method(@(x) H1_temp*x, @(x) H1_temp'*x, size(u1_temp)) ;
gamma1 = 1.9/Lips_temp ;
% ----------------------------------------------------
% Iterations
% ----------------------------------------------------
for nit1 = 1:param_dde.JU1
% erreur observation
diff = H1_temp*u1_temp - YY_temp_alpha ; %YY_temp{alpha} ;
% gradient step
grad1 =  H1_temp'*diff ;
grad2 = param_dde.nu *(u1_temp - u2_temp) ;
grad = grad1 + grad2 ;
g = u1_temp - gamma1 * grad ;
% proximity step
vr = min(max( real(g), param_dde.theta_min), param_dde.theta_max) ;
vi = min(max( imag(g), param_dde.theta_min), param_dde.theta_max) ;
u1_temp = vr + 1i*vi ;
end
% ----------------------------------------------------
dk1_temp{alpha} = transpose(u1_temp(param_dde.Ind_conj)) ; 
end
% update U1
param_dde.U1{s} = dk1_temp ;
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Update U2
% -------------------------------------------------------------------------
for s = 1:param_data.T
d1_temp = cell(param_data.na,1) ;
Supp_temp1 = param_dde.Support_Dalpha1{s} ;
dk1_temp = param_dde.U1{s} ;
parfor alpha = 1:param_data.na
d1_temp{alpha} = sparse(1,Ntot) ;
d1_temp{alpha}(Supp_temp1{alpha}) = dk1_temp{alpha}(param_dde.Ind_conj) ;
end
H_temp = cell2mat(d1_temp)*X{i} ;
YY_temp = YY2{s,i} ;
Supp_temp2 = param_dde.Support_Dalpha2{s} ;
dk2_temp = param_dde.U2{s} ;
parfor alpha = 1:param_data.na
% ----------------------------------------------------
% construction of the observation matrix
H2_temp = H_temp(:,Supp_temp2{alpha}) ;
H2_temp(alpha,:) = [] ;
% construction of data vector
YY_temp_alpha = YY_temp{alpha} ;
YY_temp_alpha(alpha,:) = [] ;
% initialization of u2
u2_temp = (dk2_temp{alpha})' ; 
% construction of sigma(u1) for regularization term
u1_temp = (dk1_temp{alpha})' ; 
% step size
Lips_temp = param_dde.nu + pow_method(@(x) H2_temp*x, @(x) H2_temp'*x, size(u1_temp)) ;
gamma2 = 1.9/Lips_temp ;
% ----------------------------------------------------
% Iterations
% ----------------------------------------------------
for nit = 1:param_dde.JU2
% gradient step
diff = H2_temp*u2_temp - YY_temp_alpha ; %YY_temp{alpha} ;
grad1 = H2_temp'*diff ;
grad2 = param_dde.nu *(u2_temp - u1_temp) ;
grad = grad1 + grad2 ;
g = u2_temp - gamma2 * grad ;
% proximity step
vr = min(max( real(g), param_dde.theta_min), param_dde.theta_max) ;
vi = min(max( imag(g), param_dde.theta_min), param_dde.theta_max) ;
u2_temp = vr + 1i*vi ;
end
% ----------------------------------------------------
dk2_temp{alpha} = (u2_temp)' ; 
end
% update U2
param_dde.U2{s} = dk2_temp ;
end
% -------------------------------------------------------------------------
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Stopping criterion
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
stop_crit_temp = cell(1,param_data.T) ;
U1 = param_dde.U1 ;
U2 = param_dde.U2 ;
parfor s = 1:param_data.T
stop_crit_temp{s} = max( ...
    norm(cell2mat(U1{s}) - cell2mat(U1_old{s}),'fro')/norm(cell2mat(U1{s}),'fro'), ...
    norm(cell2mat(U2{s}) - cell2mat(U2_old{s}),'fro')/norm(cell2mat(U2{s}),'fro') ) ;
end
stop_crit =  max( cell2mat(stop_crit_temp) ) ;
if stop_crit(end) < param_dde.tol_norm 
    disp(['D: stop crit reached, it glob = ',num2str(nit_tot)])
    disp(['D: stop crit reached, it inner = ',num2str(nbit_d)])
    break
end
% -------------------------------------------------------------------------
end



%% ************************************************************************
% *************************************************************************
% Update Image    
% *************************************************************************
% *************************************************************************

% -------------------------------------------------------------------------
% Create convolution matrix G
% -------------------------------------------------------------------------
G = create_matrix_G_D1_D2_toolbox(param_data, param_dde, param_im)  ;
B_tmp   = @(x)   G * param_im.TF(x);  
Bt_tmp  = @(x)   real(param_im.TF_adj(G' * x));
Ax = pow_method(B_tmp, Bt_tmp, size(param_im.x{1})); % norme
% -------------------------------------------------------------------------

for iter_x = 1:param_algo.Jeps
    
    for i = 1:param_im.P
    ytmp{i} = param_data.y{i}-B_tmp(xo{i}) ;
eps_old = epsilon{i} ; 
grad_eps = Bt_tmp(B_tmp(epsilon{i}) - ytmp{i}) ;
eps_tmp = epsilon{i} - (1.9/Ax) .* real(grad_eps) ; 
epsilon{i} = prox_func(eps_tmp,Ax) ;
if strcmp(param_algo.opt_dict, 'Dirac')
eps_tmp = epsilon{i}(xo{i}>0) ;
eps_tmp(eps_tmp>param_im.max_x{i}) = param_im.max_x{i}(eps_tmp>param_im.max_x{i}) ;
epsilon{i}(xo{i}>0) = eps_tmp ;
end


% -------------------------------------------------------------------------
% Stopping criterion
% -------------------------------------------------------------------------
% if iter_x>10 && ( norm(epsilon-eps_old, 'fro') / norm(epsilon, 'fro') ) < param_algo.tol_x
%     disp(['x: stop crit reached, it glob = ',num2str(nit_tot)])
%     disp(['x: stop crit reached, it inner = ',num2str(iter_x)])
%     break
% end
% -------------------------------------------------------------------------
param_im.x{i} = xo{i}+epsilon{i} ;
end
% -------------------------------------------------------------------------


%% ************************************************************************
% *************************************************************************
% Display
% *************************************************************************
% *************************************************************************
Time_tot(nit_tot+1) = toc;


data_fid =  0.5* sum( abs( B_tmp(param_im.x{i}) - param_data.y{i}).^2 ) ;
reg_x = regul_x(epsilon{i}) ; 
reg_d = 0 ;
for s = 1:param_data.T
dk1_temp = param_dde.U1{s} ;
dk2_temp = param_dde.U2{s} ;
err_temp = cell(1,param_data.na) ;
parfor alpha = 1:param_data.na
err_temp{alpha} = sum( (abs(dk1_temp{alpha} - dk2_temp{alpha})).^2 ) ;
end
reg_d = reg_d + sum( cell2mat(err_temp) ) ;
end
crit(nit_tot+1) = data_fid + param_dde.nu * reg_d ;
crit(nit_tot+1) = crit(nit_tot+1) + reg_x ;
% error_dirty(nit_tot+1) = sqrt(sum( ( reshape(Bt_tmp(B_tmp(param_im.x{i}) - param_data.y{i} ),prod(param_im.Ni-1),1)).^2 )) ;
% -------------------------------------------------------------------------
if mod(nit_tot,1)==0
disp(['Global It: ',num2str(nit_tot)])
disp(['Crit: ',num2str(crit(nit_tot+1))])
disp(['Total time: ',num2str(sum(Time_tot))])
% disp(['Error dirty image: ',num2str((error_dirty(nit_tot+1)))])
disp(['SNR dirty image: ',num2str(SNR(Bt_tmp(B_tmp(param_im.x{i})), Bt_tmp(param_data.y{i})))])
disp(['SNR rec image: ',num2str(SNR(epsilon{i},param_im.im_true{i} - param_im.xo{i}))])
disp(' ')
end
end
% -------------------------------------------------------------------------


% *************************************************************************
% *************************************************************************
% Global stopping crit.
% *************************************************************************
% *************************************************************************
if (crit(nit_tot) - crit(nit_tot+1))/crit(nit_tot+1) < param_algo.tol_crit
    disp(['Global stop crit reached, it glob = ',num2str(nit_tot)])
    break
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


end






%% Final variables

results.eps_rec = epsilon ;
results.im_rec = param_im.x ;
results.U1_rec = param_dde.U1 ;
results.U2_rec = param_dde.U2 ;



end