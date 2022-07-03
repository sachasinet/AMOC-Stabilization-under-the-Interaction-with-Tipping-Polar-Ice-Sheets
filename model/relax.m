function [relaxed]=relax(state,pars,model)
% =========================================================================
% This function is used to relax to equilibrium state, solving 0=F(x).
% It is no transient, only a newton correction.
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
% =========================================================================

    stateguess = state;  % c est le guess qui va etre corrige
    
    % =====================================================================
    % Here Newton begins
    % =====================================================================
    niters = 40;
    for j = 1:niters
        
        % We begin by computing the right hand side of F and then correct
        % it with differential terms
        rvec = -F(stateguess,pars,model);
        
        
        % -----------------------------------------------------------------
        % Tolerance check:
        % -----------------------------------------------------------------
        nrm = norm(rvec);
        if (nrm < 1e-7 * sqrt(numel(state)))
            break;
        end

        [~] = J(stateguess,pars,model);
        
        
        % and now we solve the system Jdstate=-rvec, and update the
        % unknowns
%         stateguess(end)
        
        dstateguess    =  solver(rvec,pars,model);  %(z)
        stateguess       = stateguess + dstateguess;        
        
        % now just an error message in case
        if j == niters
            % -----------------------------------------------------------------
            % Show warning (perhaps change to throwing a proper exception)
            % -----------------------------------------------------------------
            fprintf('---------NEWTON FAILED---------\n')
            fprintf(' norm(F): %1.3e, iters: %d, niters: %d \n',nrm,j,niters)
            fprintf('-------------------------------\n')
        end

        
    end
    
    relaxed = stateguess;



end