function [uab, vab, Ant] = uaub2uab_disc(na,u,v,nstep, Nx)

%-PA
% Generate u_{t,alpha,beta}, v_{t,alpha,beta} (discrete indices to get MATLAB position in the image)
% Nx:    image dimension (square image)
% nstep: number of time instants
% na:    number of antennas

Ant = cell(nstep,2) ;
uab = cell(nstep,1) ;
vab = cell(nstep,1) ;

for s = 1:nstep % [] modification possible pour éviter quelques boucles, et surtout la modification de la taille des cellules impliquées
    
    Ant{s}{1} = [] ; 
    Ant{s}{2} = [] ;
    uab{s} = [] ;
    vab{s} = [] ;
    
    for alpha = 1:na-1
        Ant{s}{1} = [Ant{s}{1} ; alpha * ones((na-alpha),1)] ; % antenna alpha
        Ant{s}{2} = [Ant{s}{2} ; (alpha+1:na)'] ;              % antenna beta
        for beta = alpha+1:na
            uab{s} = [uab{s} ; u{alpha}(s) - u{beta}(s)] ; % u_alpha - u_beta
            vab{s} = [vab{s} ; v{alpha}(s) - v{beta}(s)] ; % v_alpha - v_beta
        end
    end
    uab{s} = uab{s} + Nx/2 ; % recentering (discrete indices) to get the MATLAB indices within the image (2D indices)
    vab{s} = vab{s} + Nx/2 ; % entre 0 et N-1 -> +1 par la suite
end

% % % Uab = cell2mat(uab) ;
% % % Vab = cell2mat(vab) ;
% % % figure
% % % plot(Uab, Vab, '.')

end