function [S, dS] = sealevel(xg, pars,nondim)
% (By Erik Mulder)

    global h0 x0
    
    if nargin < 3
        nondim = false;
    end
    
    if nondim == true
        xg = xg*x0;
    end
    
    [N,ds,A,beta,a,rhog,C,V,W,m,n,g,rho_ice,rho_water, sl_A, sl_B, ...
     sl_C, sl_k, sl_Lr, sl_Li] = setvalues(pars);

    dS = -sl_A*(xg-sl_Lr).^sl_k + sl_B;
    S  = -sl_A/(sl_k+1) * (xg-sl_Lr).^(sl_k+1) + sl_B*xg + sl_C;
    
    if nondim==true
        S  = S  / h0;
        dS = dS * x0 / h0;
    end
end
