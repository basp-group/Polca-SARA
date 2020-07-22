% Script to define amplitude constraints and initialisation for DIEs/DDEs

c = floor(S2/2) + 1;
p = floor(P/2) + 1;

param_dde.theta_maxR = A_center*ones(S2, P, 4)*sqrt(F);
param_dde.theta_maxR(:,:,2) = A_off*ones(S2, P, 1)*sqrt(F);
param_dde.theta_maxR(:,:,3) = A_off*ones(S2, P, 1)*sqrt(F);
param_dde.theta_maxR(:,p,:) = 0;
param_dde.theta_maxR(c,:,2) = A*ones(1,P,1)*sqrt(F);
param_dde.theta_maxR(c,:,3) = A*ones(1,P,1)*sqrt(F);
param_dde.theta_minR = -param_dde.theta_maxR;
param_dde.theta_maxI = param_dde.theta_maxR;
param_dde.theta_minI = -param_dde.theta_maxI;

% For DIEs
for i =1:3:4
    param_dde.theta_maxR(c,p,i) = sqrt(F);
    param_dde.theta_minR(c,p,i) = sqrt(F);
    param_dde.theta_maxI(c,p,i) = 0;
    param_dde.theta_minI(c,p,i) = -param_dde.theta_maxI(c,p,i);
end

for i = 2:3
    param_dde.theta_maxR(c,p,i) = (off_diag)*sqrt(F);
    param_dde.theta_maxI(c,p,i) = 0;
    param_dde.theta_minR(c,p,i) = (off_diag)*sqrt(F);
    param_dde.theta_minI(c,p,i) = 0;
end


%%
% DDE initialisation

vec_z = zeros(S2, na, P, 1);
vec_z(c,:,:,:) = A;
p = floor(P/2) + 1;
U1 = ((randn(S2, na, P, 4) + 1i*randn(S2, na, P,4))/P)*sqrt(F);
U1r = zeros(size(U1));
U1i = U1r;

U1(:,:,:,1) = A_center*U1(:,:,:,1);
U1(:,:,:,4) = A_center*U1(:,:,:,4);
U1(:,:,:,2) = vec_z.*U1(:,:,:,2);
U1(:,:,:,3) = vec_z.*U1(:,:,:,3);

U1(:,:,p,:) = 0;
U1(c,:,p,1) = sqrt(F);
U1(c,:,p,4) = sqrt(F);
U1(c,:,p,2) = sqrt(F)*off_diag;
U1(c,:,p,3) = sqrt(F)*off_diag;


for i = 1:4
    U1r(:,:,:,i) = max(min(real(U1(:,:,:,i)), reshape(param_dde.theta_maxR(:,:,i), [S2,1,P,1])), reshape(param_dde.theta_minR(:,:,i), [S2,1,P,1]));
    U1i(:,:,:,i) = max(min(imag(U1(:,:,:,i)), reshape(param_dde.theta_maxI(:,:,i), [S2,1,P,1])), reshape(param_dde.theta_minI(:,:,i), [S2,1,P,1]));
    U1(:,:,:,i) = U1r(:,:,:,i) + 1i*U1i(:,:,:,i);
end

U2 = U1;
