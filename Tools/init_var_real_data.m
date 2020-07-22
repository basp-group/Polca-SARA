% Initialize variables

for i = 1:5

  
        
  St_im{i} = zeros(Ny,Nx);
% if i<=P
% St_im{i} = im{i};
% else
%    St_im{i} = zeros(Nx,Ny); 
% end
 
%  dual_tilde1{i} = zeros(size(param_sim_data.Psit(St_im{i})));
for k=1:Ps
 dual_var1{i,k} = zeros(size(Psit{k}(St_im{i})));
end

 dual_tilde2{i} = zeros(size(St_im{i}));
 dual_var2{i} = zeros(size(St_im{i}));
 
 
 dual_tilde3{i} = zeros(size(A{1}(St_im{i})));
 dual_var3{i} = zeros(param.M,1);
 

%  rel_sol{i} = zeros(tmax,1);
 
end
% 


if method >1
epi_h{1} = 0*[dual_var2{1}(:),dual_var2{4}(:)];
epi_h{2} = 0*[dual_var2{2}(:),dual_var2{3}(:), dual_var2{5}(:)];
end

