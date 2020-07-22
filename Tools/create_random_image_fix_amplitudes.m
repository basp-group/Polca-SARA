%% Position of sources in three levels, vary with seed

x = zeros(Nx,Ny) ;
x1 = x ;
x2 = x;

rng(seed_data)

ind = randi([c_inf c_sup], N_bright,2) ;
ind = unique(ind,'rows') ;
size_ind = max(size(ind)) ;
while size_ind<N_bright
ind_tmp = randi([c_inf c_sup], N_bright-size_ind,2) ;
ind = unique([ind ; ind_tmp],'rows') ;
size_ind = max(size(ind)) ;
end


ind_s = randi([ c_inf c_sup ], n_small_imp,2 ) ;
ind_s = unique(ind_s, 'rows') ;
[inter, ind_inter1] = intersect(ind_s, ind, 'rows') ;
if size(inter,1) >0 
ind_s(ind_inter1,:) = [] ;
end
size_ind_s = size(ind_s,1) ;

while size_ind_s < n_small_imp
ind_tmp = randi([ c_inf c_sup ], n_small_imp-size_ind_s,2 ) ;
ind_s = unique([ind_s ; ind_tmp], 'rows') ;
[inter, ind_inter1] = intersect(ind_s, ind, 'rows') ;
if size(inter,1) >0 
ind_s(ind_inter1,:) = [] ;
end
size_ind_s = size(ind_s,1) ;
end
ind_S = ind_s ;

while size_ind_s < N_small
ind_tmp = randi([ 2+min_ind , Nx-min_ind ], N_small-size_ind_s,2 ) ;
ind_s = unique([ind_s ; ind_tmp], 'rows') ;
[inter, ind_inter1] = intersect(ind_s, ind, 'rows') ;
if size(inter,1) >0 
ind_s(ind_inter1,:) = [] ;
end
size_ind_s = size(ind_s,1) ;
end
[inter, ind_inter1] = intersect(ind_s, ind_S, 'rows') ;
ind_s(ind_inter1,:) = [] ;


%%

for i = 1:N_bright
x(ind(i,1), ind(i,2)) = val(i) ;
end
H = fspecial('gaussian',[3 3],0.5);
X1 = conv2(x,H,'same') ; % Image with known bright sources
x = X1(:) ;
ind_B = x>0 ;


% Filter for fainter sources
H = fspecial('gaussian',[3 3],0.15);

for i = 1:n_small_imp
x1(ind_S(i,1), ind_S(i,2)) = val_S(i) ;
end
X2 = conv2(x1,H,'same'); % Image with unknown important sources


for i = 1:N_small-n_small_imp
x2(ind_s(i,1), ind_s(i,2)) = val_s(i) ;
end
X3 = conv2(x2,H,'same'); % Image with unkonwn faint background sources

im = X1 + X2 + X3 ;
im = (im+abs(im))./2;

% Normalised images for three levels
X1 = X1./max(im(:));
X2 = X2./max(im(:));
X3 = X3./max(im(:));


im = im./max(im(:)); % im = X1 + X2 + X3

im_app = im(:) ;
im_app(logical(1-ind_B)) = 0 ;
im_orig_app = reshape(im_app,Nx,Ny) ;
im_app = im_orig_app ;

%%
if param_im.error

im_err = randn(size(im_app));
im_err = param_im.std_err.*im_err ;

im_err = min(max(im_err,-0.1*im_app), 0.1*im_app) ;
im_app = (1+im_err).*im_app;
im_app = reshape(im_app,Nx,Ny) ;

end


figure
subplot 321, imagesc(im), axis image, colorbar
subplot 322, imagesc(log10(im)), axis image, colorbar
subplot 323, imagesc(im_app), axis image, colorbar
subplot 324, imagesc(log10(im_app)), axis image, colorbar
subplot 325, imagesc(abs(im-im_app)), axis image, colorbar
subplot 326, imagesc(log10(abs(im-im_app))), axis image, colorbar

