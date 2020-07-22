% Script to initialise the images

x0 =x_approx; % Approximated from calibration transfer + imaging
param_algo.im_approx_thresh = 0.01; % Threshold parameter


for i =1:3
    xq{i} = x_approx{i};
    x0{i}(abs(x0{i})< param_algo.im_approx_thresh*max(abs(x0{i}(:)))) = 0;
    mask_app{i} = (abs(x0{i}) > 0);
    
    param_im.min_x{i} = zeros(size(x0{i}));
    param_im.min_x{i}(mask_app{i}) = -5e-1*abs(x0{i}(mask_app{i}));
    param_im.max_x{i} = 5e-1*abs(x0{i}(abs(x0{i})>0));

    sol_mask{i} = param_im.min_x{i};
    sol_mask{i} = x0{i}(mask_app{i});
    
    epsilon{i} = zeros(size(x0{i}));
    x_th.x_approx{i} = x0{i};
    x_th.eps_true{i} = x_th.x{i}-x_th.x_approx{i};
    x_th.xo{i} = x0{i};
end
