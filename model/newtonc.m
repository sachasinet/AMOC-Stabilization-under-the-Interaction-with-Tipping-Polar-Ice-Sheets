function [state,pars,ad]=newtonc(cid,state,pars,state0,pars0,ds,model)
% Newton iteration used in the (pseudo arc length) continuation process.
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

nvar = numel(state);
neq  = numel(F(state,pars,model));
par=pars(cid);
par0=pars0(cid);

niters = 20;

ad=0;
for k=1:niters
    
    statedff = state-state0;
    pardff   = par-par0;
    
    r      = (ds*ds) - statedff'*statedff-pardff*pardff;
    
    rvec  = F(state,pars,model);
    Fx = J(state,pars,model);
    Fp = ones(1,neq)';
    
    if strcmp(model,'roothmod')
        Fp = [0;0;0;-1;0];
    end
    
    A                  =  zeros(neq+1,nvar+1);
    A(1:neq,1:nvar)    =  Fx;
    A(1:neq,nvar+1)    =  Fp;
    A(neq+1,1:nvar)  =  -2*(statedff);
    A(neq+1,nvar+1)  =  -2*(pardff);
    
    B = [-rvec;-r];
    
    nrm = norm(B);

    if (nrm < 1e-7 * sqrt(numel(state)))
        break;
    elseif nrm>1e-4* sqrt(numel(state))
        ad = 1;
    end
    
    %A = [Fx,Fp;-2*dx,-2*dp];
    %B = [-F;-r];
    
    Deltastate = A\B;
    
    state = state + Deltastate(1:nvar);
    pars(cid) = pars(cid) + Deltastate(nvar+1);
    par=pars(cid);
    
    if k == niters
        % -----------------------------------------------------------------
        % Show warning (perhaps change to throwing a proper exception)
        % -----------------------------------------------------------------
        fprintf('---------NEWTON FAILED------------------\n')
        fprintf(' norm(F): %1.3e, iters: %d, niters: %d \n',nrm,k,niters)
        fprintf('----------------------------------------\n')
    end
    
end


end
