% Create Stokes images for Cygnus A

n1 = 256;
n2 = n1;
if n1 == 256
size_ker_h = [9,9] ;
var_ker_h = 1;
else
    size_ker_h = [5,5];
    var_ker_h = 0.5;
end

h = fspecial('gaussian',size_ker_h,var_ker_h);

% x1 = fitsread('rwsol_I_15.fits');
% im{1} = imresize(x1,[n1,n2]);
% im{1}(im{1}<1e-3) = 1e-3;
% 
% x2 = fitsread('rwsol_Q_15.fits');
% im{2} = imresize(x2,[n1,n2]);
% im{2}(im{2}<1e-4) = 1e-4;
% 
% x3 = fitsread('rwsol_U_15.fits');
% im{3} = imresize(x3,[n1,n2]);
% im{3}(im{3}<1e-4) = 1e-4;
% 


% x11 = fitsread('CYG-SHI-0.65-I.fits');
% x11 = fitsread('CYG-XLO-0.20-I.fits');
x11 = fitsread('CYG-CLO-0.35-I.fits');
% x1 = x11(1229:2892,513:3840);
x1_ = x11(1229:2892,413:3940);
x1 = x11(1179:2962,283:4040);
% x1 = x11(1319:2842,513:3840);
im{1} = imresize(x1,[n1,n2]);
% im{1}(im{1}<10^(-3.2)) = 0; %1e-3; %1e-5
% im{1} = (im{1}./max(im{1}(:))).*max(x1(:));
im_ = conv2(im{1},h,'same');
im{1} = (im_./max(im_(:))).*max(x1(:));
im{1}((im{1})<10^(-3.2)) = 0;



% x21 = fitsread('CYG-SHI-0.65-Q.fits');
% x21 = fitsread('CYG-XLO-0.20-Q.fits');
 x21 = fitsread('CYG-CLO-0.35-Q.fits');
% x2 = x21(1229:2892,513:3840);
 x2 = x21(1229:2892,413:3940);
 x2 = x21(1179:2962,283:4040);
im{2} = imresize(x2,[n1,n2]);
% im{2}(abs(im{2})<10^(-3.4)) =0; %sign(im{2}(abs(im{2})<1e-5)).*1e-5; %1e-6
% im{2} = (im{2}./max(im{2}(:))).*max(x2(:));
im_ = conv2(im{2},h,'same');
im{2} = (im_./max(im_(:))).*max(x2(:));
im{2}(abs(im{2})<10^(-3.4)) = 0; %3.5
im{2}(im{1} <= 0) = 0;

% x31 = fitsread('CYG-SHI-0.65-U.fits');
% x31 = fitsread('CYG-XLO-0.20-U.fits');
x31 = fitsread('CYG-CLO-0.35-U.fits');
% x3 = x31(1229:2892,513:3840);
 x3 = x31(1229:2892,413:3940);
 x3 = x31(1179:2962,283:4040);
im{3} = imresize(x3,[n1,n2]);
% im{3}(abs(im{3})<10^(-3.35)) = 0; %sign(im{3}(abs(im{3})<1e-5)).*1e-5; %1e-6
% im{3} = (im{3}./max(im{3}(:))).*max(x3(:));
 im_ = conv2(im{3},h,'same');
im{3} = (im_./max(im_(:))).*max(x3(:));
im{3}(abs(im{3})<10^(-3.05)) = 0;
im{3}(im{1} <= 0) = 0;

%   im = create_stokes_im(im, 1, stokes_P);
im = create_stokes_cyg_a(im, stokes_P);

%%
% figure, imagesc(log10(abs(im_true{1}))), colorbar, axis image, axis off
% figure, imagesc(log10(abs(im_true{2}))), colorbar, axis image, axis off
% figure, imagesc(log10(abs(im_true{3}))), colorbar, axis image, axis off
