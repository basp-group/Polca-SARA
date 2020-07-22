% Generate image

if strcmp(param_sim_data.im_choice,'file')

im_ex = xlsread('avery_sgra_excel_file');
im = cell(3,1);

for i = 1:3
im{i} = im_ex(9:end,i+2); 

s = size(im{i},1);
im{i} = reshape(im{i},sqrt(s),sqrt(s));

if i == 1
im{i}(im{i}<0) = 0;
end

% scale image to avoid numerical errors
%      im{i} = im{i}/max(max(im{i}));
end

else
    
Nx = param_sim_data.Nx;
Ny = param_sim_data.Ny;
    

if strcmp(param_sim_data.im_type,'avery')
    ii = fitsread('avery_hm2.I.model.fits');
    qq = fitsread('avery_hm2.Q.model.fits');
    uu = fitsread('avery_hm2.U.model.fits');
    
    im{1} = ii(230:295,220:305);
    im{2} = qq(230:295,220:305);
    im{3} = uu(230:295,220:305);
    
else
    ii = fitsread('jason-j2.I.model.fits');
    qq = fitsread('jason-j2.Q.model.fits');
    uu = fitsread('jason-j2.U.model.fits');
    
    im{1} = ii(60:240,60:240);
    im{2} = qq(60:240,60:240);
    im{3} = uu(60:240,60:240);
      
end

im{1} = imresize(im{1},[Nx,Ny]);
im{2} = imresize(im{2},[Nx,Ny]);
im{3} = imresize(im{3},[Nx,Ny]);

im{1}(im{1}<0) = -im{1}(im{1}<0);

end

% P = im{2} + 1i*im{3};
% P_max = max(abs(P(:)));


