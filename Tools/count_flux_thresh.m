% Compute flux constraint error images
% P_abs = sqrt(abs(St_im{2}).^2+abs(St_im{3}).^2);

P_abs = sqrt(abs(St_im{2}).^2 + abs(St_im{3}).^2); %v12;
flux_err = P_abs - St_im{1};

mask = zeros(size(St_im{1}));

mask(flux_err>0) = flux_err(flux_err>0);
% figure,imagesc(mask)
% colorbar('FontSize',16,'FontWeight','bold')
% axis image, axis off

% if method == 1
% flux_err1_jason(k) = size(find(mask(:)),1);
% fprintf('Number of pixel with unsat. const.= %d\n', flux_err1_jason(k));
% else
%    flux_err3_jason(k) = size(find(mask(:)),1);
% fprintf('Number of pixel with unsat. const.= %d\n', flux_err3_jason(k));
% end


%% Perform thresholding
% if gen_data == 2
% dynamic_range_calc_real;
% else
%     dynamic_range_calc;
% end

g3__ = g3{1}./sigma3;
dd = norm(real(g3__),2);
res_norm(1) = dd./(sqrt(N)*evl);

mask2 = zeros(size(St_im{1}));
mask2(flux_err>(3*res_norm(1))) = flux_err(flux_err>(3*res_norm(1)));
% figure, imagesc(mask2)
% colorbar('FontSize',16,'FontWeight','bold')
% axis image, axis off


flux_thresh =  size(find(mask2(:)),1);

% if method == 1
% flux_err1_thresh_jason(k) = size(find(mask2(:)),1);
% fprintf('Number of pixel with unsat. const. after threshold= %d\n', flux_err1_thresh_jason(k));
% else
%    flux_err3_thresh_jason(k) = size(find(mask2(:)),1);
% fprintf('Number of pixel with unsat. const. after threshold= %d\n', flux_err3_thresh_jason(k));
% end 
