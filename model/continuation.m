function [curve]=continuation(cid,state,pars,parsc,k,model)
%The main function for (pseudo arc length) continuation
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

[ds,dp]=setvalues(parsc);
[M] = numel(state);
par = pars(cid);

% we begin by relaxing the initial condition to equilibrium
state = relax(state,pars,model);

curve(1,:) = [state',par];
nn=0;

fprintf('Performing continuation ')
for i=1:k
    waiting(i,k)
    
    if i == 1
        dpars=pars;
        dpars(cid)=pars(cid)+dp;
        
        Fx  = J(state,pars,model); 
        Fp = (F(state,dpars,model)-F(state,pars,model))./dp;     
        
        dx = Fx\-Fp;
        
        dx = dx./sqrt(norm(dx)^2 + 1); 
        dp = 1/sqrt(norm(dx)^2 + 1);
        
        statedot = dx;
        pardot = dp;
        
    else
        % we build the predictor
        statepast = curve(i-1,1:M)'; 
        parpast    = curve(i-1,end);
        
        statedot = (state-statepast)./dslast;
        pardot = (par-parpast)./dslast;
    end
    
    stateguess    = state+ds.*statedot;
    parguess    = par+ds.*pardot;
    
    parsguess = pars;
    parsguess(cid) = parguess;
    
    % we correct the predictor. For this we need the function, the jacobian
    %, and the derivative of the relation of the arclength
    
    statesave=state;
    parssave=pars;
    parsave = parssave(cid);
    [state,pars,ad] = newtonc(cid,stateguess,parsguess,state,pars,ds,model);
    nretry = 0;
    
    %we reduce the step if precision was not reached
    while ad == 1 && ~(i==1)
        
        if ds==parsc(1)
            nn = nn + 1;
        end

        ds = ds/10;
        stateguess    = statesave+ds.*statedot;
        parguess    = parsave+ds.*pardot;
        parsguess = pars;
        parsguess(cid) = parguess;
        [state,pars,ad] = newtonc(cid,stateguess,parsguess,statesave,parssave,ds,model);
        
        nretry=nretry+1;
        
        if nretry>5
            disp('too much retry')
            break
        end
        
    end
    
    par = pars(cid);
    curve(i+1,:)=[state',par];
    dslast = ds;
    [ds,dp]=setvalues(parsc);
    
    
end
fprintf('Continuation done \n\n')
end