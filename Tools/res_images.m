%% Plot residual images

 b_im = conv_stokes_to_bright([St_im{1}(:),St_im{2}(:),St_im{3}(:)],L,size(prev_xsol{1}));
 
 
for i = 1:4
    
     res{i} = A{i,1}(b_im{i});
     res{i} = res{i} - y(:,i);
     res_{i} = At{i,1}(res{i,1});
end

figure(600)
subplot 221, imagesc(real(res_{1}))
subplot 222, imagesc(real(res_{2}))
subplot 223, imagesc(real(res_{3}))
subplot 224, imagesc(real(res_{4}))
