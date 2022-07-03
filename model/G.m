function [out]=G(F,state,pars,statep,dt,model)
% Here I define G = F(x) - Bdx/dt used to solve transient 
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

global arr t0R t0

if strcmp(model,'S')
    %% South
    % We modify it to get G
    
    [N,ds] = setvalues(pars);
    s = (0:ds:1)';
    
    hp    =  statep(1:N);
    xgp   = statep(2*N+1);
    
    h   = state(1:N);
    xg  = state(2*N+1);
    
    %First we go for the rhs part
    out=F;
    %==========================================================================
    % Left boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    out(1) = F(1) - (h(1)-hp(1))/dt ;
    out(2) = F(2) - (h(2)-hp(2))/dt+s(2)/xg*(h(3)-h(1))/(2*ds)*(xg-xgp)/dt;
    %==========================================================================
    % Interior
    %==========================================================================
    k = 3:N-1;
    % mass conservation----------------------------------------------------
    out(k) = F(k) - (h(k)-hp(k))/dt+s(k)/xg.*(h(k+1)-h(k-1))./(2.*ds).*(xg-xgp)./dt;
    %==========================================================================
    % Right boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    out(N) =  F(N) - (h(N)-hp(N))/dt+s(N)/xg*(h(N)-h(N-1))/(2*ds)*(xg-xgp)/dt;
    
    %Now we go for the Jacobian part
    %==========================================================================
    % Left boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    arr.D1(1) = arr.D1(1)-1/dt;
    arr.D1(2) = arr.D1(2)-1/dt;
    %--------------------------------------------------------------------------
    % Ici va falloir ajouter des U(1) et des L(1) pour la seconde eq
    arr.U1(2) =  s(2)/(xg*2*ds)*((xg-xgp)/dt);
    arr.L1(2) =  -s(2)/(xg*2*ds)*((xg-xgp)/dt);
    %--------------------------------------------------------------------------
    arr.R(2)  = arr.R(2) + s(2)/(xg^2)*(h(3)-h(1))/(2*ds)*xgp/dt;
    
    %==========================================================================
    % Interior
    %==========================================================================
    k = 3:N-1;
    %----------------------------------------------------------------------
    % mass conservation
    %----------------------------------------------------------------------
    arr.D1(k)   =  arr.D1(k)-1/dt;
    %--------------------------------------------------------------------------
    % Ici va falloir ajouter des U(i) et des L(i)
    arr.U1(k) =  s(k)/(xg*2*ds)*((xg-xgp)/dt);
    arr.L1(k) =  arr.L1(k) -s(k)/(xg*2*ds)*((xg-xgp)/dt);
    %----------------------------------------------------------------------
    arr.R(k)    =  arr.R(k) + s(k)./(xg^2).*(h(k+1)-h(k-1))./(2*ds)*xgp/dt;
    
    %==========================================================================
    % arr.Right boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    arr.D1(N) = arr.D1(N)  + s(N)/(xg*2*ds)*((xg-xgp)/dt)-1/dt;
    arr.L1(N) = arr.L1(N) -s(N)/(xg*2*ds)*((xg-xgp)/dt) ;
    %--------------------------------------------------------------------------
    arr.R(N)   =  arr.R(N) + s(N)/(xg^2)*(h(N)-h(N-1))/(2*ds)*xgp/dt;
    
elseif strcmp(model,'roothmod')
    %% Rooth
    out = F - (state-statep)/dt;
    arr.JR = arr.JR-eye(5,5)./dt;
    
elseif strcmp(model,'N')
    %% North
    out = F - (state-statep)/dt;
    arr.DN = arr.DN - 1/dt;
        
elseif strcmp(model,'SRN')
    %% SRN
    [FS,parsS,FR,parsR,FN,parsN] = cut(F,pars,model);
    [stateS,~,stateR,~,stateN,~] = cut(state,pars,model);
    [statepS,~,statepR,~,statepN,~] = cut(statep,pars,model);
    
    outS = G(FS,stateS,parsS,statepS,dt,'S');
    outR = G(FR,stateR,parsR,statepR,dt,'roothmod');
    outN = G(FN,stateN,parsN,statepN,dt,'N');
    
    [NS]=setvalues(parsS);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    
    [Vol,dVol]   = VolNmod(stateN,parsN);
    [Volp,~] = VolNmod(statepN,parsN);
    
    outR(4) = outR(4) + (Vol - Volp)/dt * sigmaN;
    
    xg = stateS(2*NS+1);
    xgp = statepS(2*NS+1);
    hxg = stateS(NS);
    outR(5) = outR(5) + hxg*(xg - xgp)/dt * sigmaS*t0/t0R;
    
    [NN] = setvalues(parsN);
    
    i = 1:NN;
    arr.CRN(4,i)  = arr.CRN(4,i) + dVol(i)./dt .* sigmaN;
    
    arr.CRS(5,2*NS+1)  = hxg*1/dt*sigmaS*t0/t0R; 
    arr.CRS(5,NS)  = arr.CRS(5,NS) + (xg - xgp)/dt * sigmaS*t0/t0R;
    
    %reconstruction
    out = [outS',outR',outN']';
    
elseif strcmp(model,'SRNh')
    %% hosing
    [FS,parsS,FR,parsR,FN,parsN] = cut(F,pars,model);
    [stateS,~,stateR,~,stateN,~] = cut(state,pars,model);
    [statepS,~,statepR,~,statepN,~] = cut(statep,pars,model);
    
    outS = G(FS,stateS,parsS,statepS,dt,'S');
    outR = G(FR,stateR,parsR,statepR,dt,'roothmod');
    outN = G(FN,stateN,parsN,statepN,dt,'N');
    
    [NS]=setvalues(parsS);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    
    [Vol,dVol]   = VolNmod(stateN,parsN);
    [Volp,~] = VolNmod(statepN,parsN);
    
    outR(4) = outR(4) + (Vol - Volp)/dt * sigmaN;
    
    xg = stateS(2*NS+1);
    xgp = statepS(2*NS+1);
    hxg = stateS(NS);
    outR(5) = outR(5) + hxg*(xg - xgp)/dt * sigmaS*t0/t0R;
    
    [NN] = setvalues(parsN);
    
    i = 1:NN;
    arr.CRN(4,i)  = arr.CRN(4,i) +dVol(i)./dt .* sigmaN;
    
    arr.CRS(5,2*NS+1)  = hxg*1/dt*sigmaS*t0/t0R;
    arr.CRS(5,NS)  = arr.CRS(5,NS) + (xg - xgp)/dt * sigmaS*t0/t0R;
    
    %reconstruction
    out = [outS',outR',outN']';
    
end

end