
% n_lev : number of intensity levels
% rho : percentage for std deviation: std_lev = rho * moy_lev
% moy_lev : vector of size n_lev
%           containing the mean value of intensity sources of the 
%           different levels
% std_lev : vector of size n_lev
%           containing the std deviation for the intensity sources of the 
%           different levels  
% I_lev : vector of size n_lev
%         containing the number of sources per level
% En_lev : vector of size n_lev
%          containing the total energy of each level
% Nx : size image
% c_inf, c_sup : sub part of the image where bright sources belong


%% Testing parameters


% % % close all
% % % clear all
% % % clc
% 
% ws = 0 ;
% seed_data = 0 ;
% Nx = 128; Ny = 128 ;
% 
% gen_amp_sources = 1 ;
% n_lev = 3 ;
% I_lev = [10, 10, 200] ;
% En_lev = [10, 1, 1e-6] ;
% rho = 0.1 ;
% 
% size_ker_h = [3,3] ;
% var_ker_h = 0.5 ;
% h = fspecial('gaussian',size_ker_h,var_ker_h);



%%


rng(seed_data)


%% Generate positions


c = floor(Nx/2)+1 ;

c_inf = c- floor(Nx/(3)) ;
c_sup = c+ floor(Nx/(3)) ;

min_ind = 3 ;

ind_tot = [] ;

for l = 1:n_lev

if l <=2 
window_sources = [c_inf c_sup] ;
else
window_sources = [ 2+min_ind , Nx-min_ind ] ;
end

ind = randi(window_sources, I_lev(l),2) ;
ind = unique(ind,'rows') ;
if l>1
[inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
if size(inter,1) >0 
ind(ind_inter,:) = [] ;
end
end
size_ind = max(size(ind)) ;
while size_ind<I_lev(l)
ind_tmp = randi(window_sources, I_lev(l)-size_ind,2) ;
ind = unique([ind ; ind_tmp],'rows') ;
[inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
if size(inter,1) >0 
ind(ind_inter,:) = [] ;
end
size_ind = max(size(ind)) ;
end


if l == 1
ind_lev{l} = ind ;
x_tmp = zeros(Nx, Ny) ;
for i =1:I_lev(l)
x_tmp(ind(i,1), ind(i,2)) = 1 ;
end
size_ker_tmp = size_ker_h + ( floor(size_ker_h./2) .* 2 ) ;
h_tmp = fspecial('gaussian',size_ker_h+2,var_ker_h);
x_tmp = conv2(x_tmp,h_tmp,'same') ;
x_tmp = x_tmp(:) ;
ind_tmp = find(x_tmp>0) ;
ind_tot = convert_1D_2D( ind_tmp, Nx) + (Nx/2+1) ;
ind_tot = [ind_tot(:,2), ind_tot(:,1)] ;
else
ind_lev{l} = ind ;
ind_tot = [ind_tot ; ind] ;
end
    

end






%% Generate amplitudes
if gen_amp_sources == 1
moy_lev = zeros(n_lev,1) ;
std_lev = zeros(n_lev,1) ;
min_lev = zeros(n_lev,1) ;
max_lev = zeros(n_lev,1) ;
amp_lev = cell(1,n_lev) ;
for l = 1:n_lev

moy_lev(l) = En_lev(l) / sqrt( I_lev(l) * (rho^2+1) ) ;
std_lev(l) = rho * moy_lev(l) ;

amp_lev{l}(1:I_lev(l)) = std_lev(l) * randn(I_lev(l),1) + moy_lev(l) ;

min_lev(l) = min(amp_lev{l}) ;
max_lev(l) = max(amp_lev{l}) ;

end
end



%% Create images


im = zeros(Nx,Ny) ;
ind_nz = cell(1,n_lev) ; % storing indices of non-zero elements
x_level = cell(length(n_lev),1) ;

for l = 1:n_lev

x_level{l} = zeros(Nx,Ny) ;
for i = 1:I_lev(l)
x_level{l}(ind_lev{l}(i,1), ind_lev{l}(i,2)) = amp_lev{l}(i) ;
end
x_level{l} = conv2(x_level{l},param_im.h,'same') ;
energy_temp = sqrt(sum(x_level{l}(:).^2)) ;
alpha = En_lev(l) / energy_temp ;
x_level{l} = alpha*x_level{l} ;

disp(['level: ',num2str(l)])
disp(['Theoretical energy: ', num2str(En_lev(l))])
disp(['True energy: ', num2str( sqrt(sum( x_level{l}(:).^2 )) )])
disp(['number of sources: ',num2str(I_lev(l))])
disp(' ')
if ws == 0
figure, imagesc(x_level{l}), axis image, colorbar
end


im = im+x_level{l} ;
end
param_im.x_level = x_level ;
param_im.ind_lev = ind_lev ;

im_app = x_level{1} ;

if param_im.error % introduce error in bright sources
im_err = randn(size(im_app)) ;
im_err = param_im.std_err.*im_err ;

im_err = min(max(im_err,-0.1*im_app), 0.1*im_app) ;
im_app = (1+im_err).*im_app;
im_app = reshape(im_app,Nx,Ny) ;
else
im_app = im_app ;
end


% 
% figure
% subplot 221, imagesc(log10(x_level{1})), axis image, colorbar
% subplot 222, imagesc(log10(x_level{2})), axis image, colorbar
% subplot 223, imagesc(log10(x_level{1}+x_level{2})), axis image, colorbar


if ws == 0
figure
subplot 321, imagesc(im), axis image, colorbar
subplot 322, imagesc(log10(im)), axis image, colorbar
subplot 323, imagesc(im_app), axis image, colorbar
subplot 324, imagesc(log10(im_app)), axis image, colorbar
subplot 325, imagesc(abs(im-im_app)), axis image, colorbar
subplot 326, imagesc(log10(abs(im-im_app))), axis image, colorbar
end



