% Dynamic range

% load('/Users/jasleenbirdi/Documents/PhD/Polarization_imaging/Results/meas_operator.mat')

n1 = size(im{1},1);
n2 = size(im{1},2);
N = n1*n2;

res = cell(P,1);

p=1;

         x = St_im{p};
% x= im{p};

%     x = St_im{p};
    % x = im{P};
    %     num = sqrt(N)* evl.^2 * max(x(:));
    
    peak(p) = max(x(:));
    rms_n(p) = rms(St_im{p}(:));
    num(p) = sqrt(N)* evl * max(x(:));
    
    den_term = real(At(Gw'*((Gw*A(x)-y{1,n_test}{p,1}))));
    res{p} = den_term;
    den(p) = norm(den_term,2);
    res_norm(p) = den(p)./(sqrt(N)* evl);
    DR(p) = num(p)./den(p);
    
    DR_peak_rms(p) = peak(p)./rms_n(p);




