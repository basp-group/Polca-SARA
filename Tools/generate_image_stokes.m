function [x0, eps_true, im_min, im_max, min_xo,  E1, E2, tau] = generate_image_stokes(kappa, im)

SNR =@(x, xtrue) 20 * log10( sqrt( sum( xtrue(:).^2 ) / sum( (x(:)-xtrue(:)).^2 ) ) );

%% Load image
% 
% im = fitsread([im_choice, '.fits']);
% if any(N < size(im))
%     im = imresize(im, N);
% end
% im = (im+abs(im))./2;
% im = im./max(im(:));

%% Generate approximated image from the true image

% energy of the images image
E_im = sqrt(sum(im(:).^2));

E_im_app = sqrt(1/(1+kappa^2))*E_im;
E_im_b = kappa*E_im_app;

% generate approximation
tau_min = 0; tau_max = max(im(:)); %1;
E_max = E_im; E_min = 0;
diff_sup = E_im_app - E_max;
diff_inf = E_im_app - E_min;
if diff_sup > diff_inf
    tau = tau_max; x0 = zeros(size(im));
else
    tau = tau_min; x0 = im;
end

cond_stop = 1;
k = 0;
k_max = 100;
tol = 1e-5;

disp(['Find approximated image with theoretical energy = ',num2str(E_im_app)])
disp(['Global energy of the image : ',num2str(E_im)])
disp(['Theoretical energy of the background image : ',num2str(E_im_b)])
disp(' ')

while cond_stop>tol && k<k_max
    k = k+1;
    
    tau_tmp = (tau_min + tau_max)/2;
    
    im_app_tmp = im;
    im_app_tmp(abs(im_app_tmp)<tau_tmp) = 0;
    
    E_tmp = sqrt(sum(im_app_tmp(:).^2));   
    
    if E_im_app > E_tmp
        tau_max = tau_tmp;
        E_min = E_tmp;
        cond = E_im_app - E_min;
        diff_inf = cond;
        if cond < diff_sup
            x0 = im_app_tmp;
            tau = tau_max;
        end
    else
        tau_min = tau_tmp;
        E_max = E_tmp;
        cond = E_max - E_im_app;
        diff_sup = cond;
        if cond < diff_inf
            x0 = im_app_tmp;
            tau = tau_min;
        end
    end
    
    cond_stop = abs(E_im_app - E_tmp);
    
end

if (norm(x0-im) == 0) || (norm(x0) == 0)
    x0 = im_app_tmp;
end
    

eps_true = im - x0;

E1 = sqrt( sum(x0(:).^2) );
E2 = sqrt( sum(eps_true(:).^2) );
disp('********************************************')
disp(['Threshold parameter tau = ',num2str(tau)])
disp(['Energy of the approximated image = ',num2str(sqrt( sum(x0(:).^2) ))])
disp(['Energy of the background image = ',num2str(sqrt( sum(eps_true(:).^2) ))])
disp(['SNR approximated image = ',num2str(SNR(x0,im))])
disp('********************************************')

%% Return useful parameters
min_x0 = min(x0(:)>0);
im_min = 0;
im_max = Inf;
mask_x0 = (x0>0);
min_xo = min(x0(:)>0) ;

end