function [mu,save_mu,save_crit] = golden_search(crit,mu_inf,mu_sup)


% mu_inf = 0.1;
% mu_sup = 10;
tau     = 0.38197;
tol     = 1e-3;
maxiter = -2.078*log(tol);

mu_1 = (1-tau)*mu_inf + tau*mu_sup;
mu_2 = tau*mu_inf + (1-tau)*mu_sup;

%disp('Initialisation');
disp(['it = ',num2str(0),' | lambda = ',num2str(mu_1)]);
disp(['it = ',num2str(0),' | lambda = ',num2str(mu_2)]);

 
J_1 = feval(crit,mu_1);
J_2 = feval(crit,mu_2);

save_mu(1) = mu_1;
save_mu(2) = mu_2;
save_crit(1) = J_1;
save_crit(2) = J_2;

iter = 1;
while iter<maxiter
    if J_1>=J_2
        mu_inf = mu_1;  mu_1   = mu_2;  J_1 = J_2;
        mu_2 = tau * mu_inf + (1-tau)*mu_sup;
        disp(['it = ',num2str(iter),' | lambda = ',num2str(mu_2)]);
        J_2 = feval(crit,mu_2);
        save_mu = [save_mu(:);mu_2];
        save_crit = [save_crit(:);J_2];

        % Jsave(iter) = J_2;
    else
        mu_sup = mu_2;
        mu_2   = mu_1;
        J_2 = J_1;
        mu_1  = (1-tau)*mu_inf + tau*mu_sup;
        disp(['it = ',num2str(iter),' | lambda = ',num2str(mu_1)]);
        J_1 = feval(crit,mu_1);
        save_mu = [save_mu(:);mu_1];
        save_crit = [save_crit(:);J_1];
       % Jsave(iter) = J_1;
    end
   
    iter = iter+1;
end

if J_1>J_2,
    mu = mu_2;J = J_2;
else
    J = J_1; mu = mu_1;
end
% figure
% plot(-Jsave); 