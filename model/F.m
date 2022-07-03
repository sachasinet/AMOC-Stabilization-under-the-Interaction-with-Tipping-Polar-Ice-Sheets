function [out] = F(sol,pars,model)
%==========================================================================
% Create right hand side of the discretization.
%
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
%==========================================================================
global t0 t0R

if strcmp(model,'S')
    % Modified by Sacha Sinet, from Eric Mulder
    %% South
    [N,ds,A,beta,a,rhog,C,V,W,m,n,g,rho_ice,rho_water,~, ~, ...
        ~, ~, ~, ~, Vb, Wb, ~, ~, ~, ~, ~, btrs] = setvalues(pars);
    
    A = exp(-beta)*A;
    s = (0:ds:1)';
    a = a*ones(N,1);
    %==========================================================================
    % STATIONARY
    %==========================================================================
    out  = -ones(2*N+1,1);
    h    =  sol(1:N);
    u    =  sol(N+1:2*N);
    xg   =  sol(2*N+1);
    
    %==========================================================================
    % Bedrock, non-dimensionalized
    %==========================================================================
    [b]   =  bed(s,xg,true);
    %==========================================================================
    % Sea level, non-dimensionalized
    %==========================================================================
    [SL, ~] = sealevel(xg,pars,true);
    %==========================================================================
    % Left boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    out(1) = -h(1)*(u(1)+u(2))/(2*xg*ds)       + a(1);
    out(2) = -(1/xg)*h(2)*(u(2) + u(1))/(2*ds) + a(2);
    
    % momentum conservation------------------------------------------------
    ud1 = u(2)-u(1);
    ud2 = u(1)+u(1);
    out(N+1) = 2*V*A^(-1/n)/(xg*ds)^(1+1/n) * (...
        h(2)*abs(ud1)^(1/n-1)*ud1-h(1)*abs(ud2)^(1/n-1)*ud2)...
        - C*abs(u(1))^(m-1)*u(1) ...
        - W*rhog*(h(1)+h(2))*(1/(2*xg*ds))*(h(2)-b(2)-h(1)+b(1));
    %==========================================================================
    % Interior
    %==========================================================================
    i = 3:N-1;
    % mass conservation----------------------------------------------------
    um1    = (u(i)  +u(i-1))/2;
    um2    = (u(i-1)+u(i-2))/2;
    out(i) = -(1/xg).*(h(i).*(um1)-h(i-1).*(um2))/(ds) + a(i);
    i = 2:N-1;
    % momentum conservation------------------------------------------------
    ud1 = u(i+1)-u(i);
    ud2 = u(i)-u(i-1);
    out(i+N) = 2.*V.*A.^(-1/n)./(xg.*ds).^(1+1/n) .* (...
        h(i+1).*abs(ud1).^(1/n-1).*ud1-h(i).*abs(ud2).^(1/n-1).*ud2)...
        - C.*abs(u(i)).^(m-1).*u(i) ...
        - W*rhog.*(h(i)+h(i+1)).*(1./(2.*xg.*ds)).*(h(i+1)-b(i+1)-h(i)+b(i));
    %==========================================================================
    % Right boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    um1    =  (u(N)  +u(N-1))/2;
    um2    =  (u(N-1)+u(N-2))/2;
    out(N) =  -(1/(xg))*(h(N)*(um1)-h(N-1)*(um2))/(ds) + a(N);
    
    % momentum conservation----------------------------------------------------
    out(2*N) = -2*Vb*A^(-1/n)/(xg*ds)^(1/n)*abs(u(N)-u(N-1))^(1/n-1)*(u(N)-u(N-1))...
        +Wb*(btrs/2)*(1-rho_ice/rho_water)*rho_water*g*(b(N) + SL);
    
    %==========================================================================
    % Auxiliary
    %==========================================================================
    % extrapolation of the flotation criterium
    out(2*N+1) = 3*h(N)-h(N-1) - 2*(rho_water/rho_ice)*(b(N) + SL);
    
elseif strcmp(model,'roothmod')
    %% Rooth
    out=zeros(5,1);
    
    [~,tau1R,tau2R,tau3R,F1R,F3R]=setvalues(pars);
    
    T1=sol(1);
    T2=sol(2);
    T3=sol(3);
    
    ad=advectionm(sol,pars);
    
    out(1) = ad(1) + tau1R - T1;
    out(2) = ad(2) + tau2R - T2;
    out(3) = ad(3) + tau3R - T3;
    out(4) = ad(4) - F1R;
    out(5) = ad(5) - F3R;
    
elseif strcmp(model,'N')
    %% North radially symmetric
    [NN,drN,A0,n] = setvalues(pars);
    
    h    =  sol(1:NN);
    hend = 0;
    
    rN   =  -1:drN:1;
    rN = rN(length(rN)/2+1:end-1);
    rNend = 1;
    
    
    [b,bend] = bedN(rN,pars);
    H = b+h;
    Hend = bend; %To be changed if it is not true anymore, for non flat bed
    
    Deff = zeros(NN,1);
    Fluxeff  = zeros(NN,1);
    dFlux = zeros(NN,1);
    out  =  ones(NN,1);
    
    for i = 1:NN-1
        Deff(i) = A0*(rN(i)+rN(i+1))/2*((h(i)+h(i+1))/2)^(n+2)*(abs(H(i+1)-H(i))/(drN))^(n-1);
    end
    Deff(NN) = A0*(rN(NN)+rNend)/2*((h(NN)+hend)/2)^(n+2)*(abs(Hend-H(NN))/(drN))^(n-1);
    
    
    Fluxeff0 = 0;
    for i = 1:NN-1
        Fluxeff(i) = -Deff(i)*(H(i+1)-H(i))/(drN);
    end
    Fluxeff(NN) = -Deff(NN)*(Hend-H(NN))/(drN);
    
    dFlux(1) = 1/(rN(1)*drN)*(Fluxeff(1)-Fluxeff0);
    for i = 2:NN
        dFlux(i) = 1/(rN(i)*drN)*(Fluxeff(i)-Fluxeff(i-1));
    end
    
    for i = 1:NN
        out(i) = -dFlux(i)+accNHA(sol,pars,i);
    end
    
elseif strcmp(model,'SRN')
    %% SRN
    [solS,parsS,solR,parsR,solN,parsN] = cut(sol,pars,model);
    
    [NS] = setvalues(parsS);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    [~,~,~,~,aN,~,~,dadtN,~,DeltaTN] = setvalues(parsN);
    
    outR = F(solR,parsR,'roothmod');
    
    % the coupling
    uxg = solS(2*NS);
    hxg = solS(NS);
    
    TS  = solR(3);
    parsS(3) = Am(TS);
    
    
    outS = F(solS,parsS,'S');
    outN = F(solN,parsN,'N');
    
    outR(4) = outR(4) - sigmaN * pi* (aN+dadtN*DeltaTN)*lenN(solN,parsN)^2;
    outR(5) = outR(5) - hxg*sigmaS*uxg;
    
    %the timescaling
    t0mis = t0;
    t0Ro   = t0R;
    
    outS(1:NS) = outS(1:NS)*t0Ro/t0mis;
    
    %reconstruction
    out = [outS',outR',outN']';
    
elseif strcmp(model,'SRNh')
    %% For hosing
    [solS,parsS,solR,parsR,solN,parsN,newpar] = cut(sol,pars,model);
    
    [NS] = setvalues(parsS);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    [~,~,~,~,aN,~,~,dadtN,~,DeltaTN] = setvalues(parsN);
    
    outR = F(solR,parsR,'roothmod');
    
    % the coupling
    uxg = solS(2*NS);
    hxg = solS(NS);
    
    TS  = newpar;
    parsS(3) = Am(TS);
    
    
    outS = F(solS,parsS,'S');
    outN = F(solN,parsN,'N');
    
    outR(4) = outR(4) - sigmaN * pi* (aN+dadtN*DeltaTN)*lenN(solN,parsN)^2;
    outR(5) = outR(5) - hxg*sigmaS*uxg;
    
    %the timescaling
    t0mis = t0;
    t0Ro   = t0R;
    
    outS(1:NS) = outS(1:NS)*t0Ro/t0mis;
    
    %reconstruction
    out = [outS',outR',outN']';
    
    
end



end
