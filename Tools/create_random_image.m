function [param_im] = create_random_image(param_im)

% for the convolution to have Gaussian sources instead of just point
% sources
size_ker_h = [3,3] ; 
var_ker_h = 0.5 ;
h = fspecial('gaussian',size_ker_h,var_ker_h);


%% Generate positions


c = floor(param_im.Ni(1)/2)+1 ;
c_inf = c- floor(param_im.Ni(1)/(3)) ;
c_sup = c+ floor(param_im.Ni(1)/(3)) ;
min_ind = 3 ;
ind_tot = [] ;

for l = 1:param_im.n_lev

if l <=2 
window_sources = [c_inf c_sup] ;
else
window_sources = [ 2+min_ind , param_im.Ni(1)-min_ind ] ;
end

ind = randi(window_sources, param_im.I_lev(l),2) ;
ind = unique(ind,'rows') ;

if l>1
[inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
if size(inter,1) >0 
ind(ind_inter,:) = [] ;
end
end
size_ind = max(size(ind)) ;

while size_ind<param_im.I_lev(l)
ind_tmp = randi(window_sources, param_im.I_lev(l)-size_ind,2) ;
ind = unique([ind ; ind_tmp],'rows') ;
if l > 1
[inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
if size(inter,1) >0 
ind(ind_inter,:) = [] ;
end
end
size_ind = max(size(ind)) ;
end


if l == 1
ind_lev{l} = ind ;
x_tmp = zeros(param_im.Ni) ;
for i =1:size(ind,1) %param_im.I_lev(l)
x_tmp(ind(i,1), ind(i,2)) = 1 ;
end
h_tmp = fspecial('gaussian',size_ker_h+2,var_ker_h);
x_tmp = conv2(x_tmp,h_tmp,'same') ;
x_tmp = x_tmp(:) ;
ind_tmp = find(x_tmp>0) ;
ind_tot = convert_1D_2D( ind_tmp, param_im.Ni(1)) + (param_im.Ni(1)/2+1) ;
ind_tot = [ind_tot(:,2), ind_tot(:,1)] ;
else
ind_lev{l} = ind ;
ind_tot = [ind_tot ; ind] ;
end


% 
% while size_ind<param_im.I_lev(l)
% ind_tmp = randi(window_sources, param_im.I_lev(l)-size_ind,2) ;
% ind = unique([ind ; ind_tmp],'rows') ;
% [inter, ind_inter] = intersect(ind, ind_tot, 'rows') ;
% if size(inter,1) >0 
% ind(ind_inter,:) = [] ;
% end
% size_ind = max(size(ind)) ;
% end

% 
% if l == 1
% ind_lev{l} = ind ;
% x_tmp = zeros(param_im.Ni) ;
% for i =1:param_im.I_lev(l)
% x_tmp(ind(i,1), ind(i,2)) = 1 ;
% end
% h_tmp = fspecial('gaussian',size_ker_h+2,var_ker_h);
% x_tmp = conv2(x_tmp,h_tmp,'same') ;
% x_tmp = x_tmp(:) ;
% ind_tmp = find(x_tmp>0) ;
% ind_tot = convert_1D_2D( ind_tmp, param_im.Ni(1)) + (param_im.Ni(1)/2+1) ;
% ind_tot = [ind_tot(:,2), ind_tot(:,1)] ;
% else
% ind_lev{l} = ind ;
% ind_tot = [ind_tot ; ind] ;
% end
    

end



%% Generate amplitudes


mean_lev = zeros(param_im.n_lev,1) ; 
std_lev = zeros(param_im.n_lev,1) ;
amp_lev = cell(1,param_im.n_lev) ;

for l = 1:param_im.n_lev
mean_lev(l) = param_im.En_lev(l) / sqrt( param_im.I_lev(l) * (param_im.rho^2+1) ) ;
std_lev(l) = param_im.rho * mean_lev(l) ;

amp_lev{l}(1:param_im.I_lev(l)) = std_lev(l) * randn(param_im.I_lev(l),1) + mean_lev(l) ;
end



%% Create images

im{1} = zeros(param_im.Ni) ;
x_level = cell(length(param_im.n_lev),param_im.P) ;

for l = 1:param_im.n_lev

x_level{l,1} = zeros(param_im.Ni) ;
for i = 1:param_im.I_lev(l)
x_level{l,1}(ind_lev{l}(i,1), ind_lev{l}(i,2)) = amp_lev{l}(i) ;
end
x_level{l,1} = conv2(x_level{l,1},h,'same') ;
energy_temp = sqrt(sum(x_level{l,1}(:).^2)) ;
alpha = param_im.En_lev(l) / energy_temp ;
x_level{l,1} = alpha*x_level{l,1} ;
im{1} = im{1}+x_level{l,1} ;
end

% save the different levels
param_im.x_level = x_level ;
% save the source positions
param_im.ind_lev = ind_lev ;

% define image containing bright sources
xo{1} = x_level{1,1} ;
% 
% if mod(param_im.Ni(1),2)==0 % to have an odd image
% im{1} = im{1}(1:end-1,1:end-1) ;
% xo{1} = xo{1}(1:end-1,1:end-1) ;
% end


%% Added by jb -- Stokes images
%%%%%%%

im{2} = zeros(param_im.Ni);
im{3} = zeros(param_im.Ni);

% Number of polarized sources in each level
 Np = param_im.Np;

for l = 1:param_im.n_lev
x_level{l,2} = zeros(param_im.Ni) ;
x_level{l,3} = zeros(param_im.Ni) ;

rand_ind{l} = datasample(ind_lev{l},Np(l),'Replace',false);

for i = 1:Np(l)

% st_im{2} = im.*(2.*rand(1,1) -1); % Numbers within the range [-1,1]
% st_im{3} = im.*(2.*rand(1,1) -1);

amp2_lev{l}(i) = amp_lev{l}(i).*(3.*rand(1,1) -1) ;
amp3_lev{l}(i) = amp_lev{l}(i).*(3.*rand(1,1) -1) ;


if sqrt(amp2_lev{l}(i).^2 + amp3_lev{l}(i).^2) > amp_lev{l}(i)
    i = i-1;
else
    x_level{l,2}(rand_ind{l}(i,1), rand_ind{l}(i,2)) = amp2_lev{l}(i) ;
    x_level{l,3}(rand_ind{l}(i,1), rand_ind{l}(i,2)) = amp3_lev{l}(i) ;

end
end
x_level{l,2} = conv2(x_level{l,2},h,'same') ;
x_level{l,3} = conv2(x_level{l,3},h,'same') ;
im{2} = im{2} + x_level{l,2} ;
im{3} = im{3} + x_level{l,3} ;
end

param_im.x_level = x_level ;
% save the source positions
param_im.ind_lev = ind_lev ;

% define image containing bright sources
xo{2} = x_level{1,2};
xo{3} = x_level{1,3};

% if mod(param_im.Ni(1),2)==0 % to have an odd image
% im{2} = im{2}(1:end-1,1:end-1) ;
% xo{2} = xo{2}(1:end-1,1:end-1) ;
% 
% im{3} = im{3}(1:end-1,1:end-1) ;
% xo{3} = xo{3}(1:end-1,1:end-1) ;
% end

%%%%%%%
%%
param_im.im_true = im ;
param_im.x = im ;
param_im.xo = xo ;
param_im.min_xo{1} = min(param_im.xo{1}(:)>0) ;
param_im.min_xo{2} = min(param_im.xo{2}(:)>0) ;
param_im.min_xo{3} = min(param_im.xo{3}(:)>0) ;

param_im.eps_true{1} = param_im.im_true{1}-param_im.xo{1} ;
param_im.eps_true{2} = param_im.im_true{2}-param_im.xo{2} ;
param_im.eps_true{3} = param_im.im_true{3}-param_im.xo{3} ;


param_im.min = 0 ;
param_im.max = Inf ;
% param_im.mask_xo = (xo>0) ;

end
