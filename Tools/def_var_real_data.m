% Define variables for the algorithm
Ps = param_pd.Ps;
St_im = cell(P+2,1); %Stokes images: I,Q,U and constraint variables Z1, Z2

% intermediate dual variables
dual_tilde1 = cell(P,Ps); %for sparsity
dual_tilde2 = cell(P+2,1); 
dual_tilde3 = cell(4,R); %for data fidelity
dual_tilde4 = cell(P,1);

% dual variables
dual_var1 = cell(P,Ps); %for sparsity
dual_var2 = cell(P+2,1); 
dual_var3 = cell(4,R); %for data fidelity
dual_var4 = cell(P,1);


% To store the norms
norm1 = cell(P,Ps); %for sparsity
dummy1 = cell(P+1,1);
dummy2 = cell(P+1,1);
dummy3 = cell(P+1,1);
dummy4 = cell(P+1,1);
St_im_old = cell(P+2,1);
epi_h = cell(2,1);
N = param.Nx*param.Ny;

%% Initialize variables

for i = 1:3
    St_im{i} = epsilon{i};
        %  dual_tilde1{i} = zeros(size(param_sim_data.Psit(St_im{i})));
    for k=1:Ps
        Psitw{k} = param_l1.Psit{k};
            Psiw{k} = param_l1.Psi{k};
        dual_var1{i,k} = zeros(size(Psitw{k}(St_im{i})));
    end
    
    dual_tilde2{i} = zeros(size(St_im{i}));
    dual_var2{i} = zeros(size(St_im{i}));
end

for i =4:5
    St_im{i} = zeros(param.Ny,param.Nx);
    dual_tilde2{i} = zeros(size(St_im{i}));
    dual_var2{i} = zeros(size(St_im{i}));
end


param_dual.dual_var1 = dual_var1;
param_dual.dual_var2 = dual_var2;
param_dual.dual_tilde2 = dual_tilde2;
param_dual.Psitw = Psitw;
param_dual.Psiw = Psiw;
param_dual.St_im = St_im;

