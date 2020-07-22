
function [dual_var1, g1, norm_l1] = run_par_waverec2(dual_var1,Psit, Psi, prev_xsol, sigma1, st_im, eta,psit_xo)

dual_tilde1 = dual_var1 + Psit(prev_xsol);
T = eta/sigma1;
soft_th = sign(dual_tilde1+psit_xo).*max(abs(dual_tilde1+psit_xo)-T,0); % Modified to incorporate xo (known bright sources image)
dual_var1 = dual_tilde1-(soft_th-psit_xo);
% soft_th = sign(dual_tilde1).*max(abs(dual_tilde1)-T,0); % Modified to incorporate xo (known bright sources image)
% dual_var1 = dual_tilde1-soft_th;

g1 = (Psi(dual_var1));

% local L1 norm of current solution
norm_l1 = sum(abs(Psit(st_im)));
end
