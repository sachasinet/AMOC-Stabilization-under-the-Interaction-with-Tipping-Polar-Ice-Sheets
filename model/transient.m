function [trans,normG]=transient(state,pars,dt,theta,model,it)
%==========================================================================
% We perform a transient from the dynamical equations Bdx/dt = F(x) 
% (implicit in time)
% 
% Author: Sacha Sinet, 2014-2017, contact -> s.a.m.sinet@uu.nl
%==========================================================================

if nargin < 6
    it=0; 
end

% =========================================================================
% Now we begin the transient computation
% =========================================================================
%we set the initial states

statep     = state; 
stateguess = statep; 

rhsp = F(statep,pars,model);
% And Newton begins
niters = 40;
for j = 1:niters

    rhs  = F(stateguess,pars,model);
    
    [~] = J(stateguess,pars,model);

    rvec = Euler(rhs,rhsp,theta); 
    
    [rvec] = G(rvec,stateguess,pars,statep,dt,model);
    
    % -----------------------------------------------------------------
    % Tolerance check:
    % -----------------------------------------------------------------
    nrm = norm(rvec);
    if (nrm < 1e-7 * sqrt(numel(state)))
        probe(stateguess,statep,pars,dt,model,it); 
        break;
    end

    rvec = -rvec;
    
    dstateguess    =  solver(rvec,pars,model);
    
    stateguess     = stateguess + dstateguess;
    
    % now just an error message in case
    if j == niters
        % -----------------------------------------------------------------
        % Show warning (perhaps change to throwing a proper exception)
        % -----------------------------------------------------------------
        fprintf('---------NEWTON FAILED---------\n')
        fprintf(' norm(G): %1.3e, iters: %d, niters: %d \n',nrm,j,niters)
        fprintf('-------------------------------\n')
    end
end

trans = stateguess;
normG = nrm;



end