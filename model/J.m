function [Jac] = J(sol,pars,model)
%==========================================================================
% Calculate the jacobian for the state given by sol and pars
%
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
%==========================================================================
global arr t0 t0R
Jac=0;

if strcmp(model,'S')
    %Modified by Sacha Sinet, from Eric Mulder
    %% South
    [N,ds,A,beta,~,rhog,C,V,W,m,n,g,rho_ice,rho_water,~, ~, ...
        ~, ~, ~, ~, Vb, Wb, ~, ~, ~, ~, ~, btrs] = setvalues(pars);
    
    A = exp(-beta)*A;
    s = (0:ds:1)';
    %==========================================================================
    % STATIONARY
    %==========================================================================
    h   = sol(1:N);
    u   = sol(N+1:2*N);
    xg  = sol(2*N+1);
    
    %==========================================================================
    % Bedrock, non-dimensionalized
    %==========================================================================
    [b,db]   =  bed(s,xg,true);
    %==========================================================================
    % Sea level, non-dimensionalized
    %==========================================================================
    [~, dSL] = sealevel(xg, pars, true);
    %==========================================================================
    % Left boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    arr.D1(1) = -(u(1)+u(2))/(2*xg*ds);
    arr.D1(2) = -(1/xg)*(u(2)+u(1))/(2*ds);
    %--------------------------------------------------------------------------
    % two different speeds appearing there, for i=1 and 2
    arr.D2(1) = -h(1)*(1)/(2*xg*ds);
    arr.D2(2) = -(1/xg)*h(2)*(1)/(2*ds);
    arr.U2(1) = -h(1)*(1)/(2*xg*ds);
    arr.L2(2) = -(1/xg)*h(2)*(1)/(2*ds);
    %--------------------------------------------------------------------------
    arr.R(1)  = h(1)*(u(1)+u(2))/(2*xg^2*ds);
    arr.R(2)  = (1/xg^2)*h(2)*(u(2) + u(1))/(2*ds);
    % momentum conservation----------------------------------------------------
    ud1       =   u(2)-u(1);
    ud2       =   u(1)+u(1);
    arr.D3(1) =   2*V*A^(-1/n)/(xg*ds)^(1+1/n) * (-abs(ud2)^(1/n-1)*ud2)...
        -W*rhog/(2*xg*ds)*(h(2)-b(2)-h(1)+b(1))...
        +W*rhog*(h(1)+h(2))*(1/(2*xg*ds));
    arr.U3(1) =   2*V*A^(-1/n)/(xg*ds)^(1+1/n)*(abs(ud1)^(1/n-1)*ud1)...
        - W*rhog/(2*xg*ds)*(h(2)-b(2)-h(1)+b(1))...
        - W*rhog*(h(1)+h(2))*(1/(2*xg*ds));
    %--------------------------------------------------------------------------
    arr.D4(1) =  2*V*A^(-1/n)/(xg*ds)^(1+1/n) * (...
        -h(2)*((1/n-1)*abs(ud1)^(1/n-2)*ud1+abs(ud1)^(1/n-1))...
        -h(1)*(2*(1/n-1)*abs(ud2)^(1/n-2)*ud2+2*abs(ud2)^(1/n-1)))...
        -(m-1)*C*abs(u(1))^(m-2)*u(1)-C*abs(u(1))^(m-1);
    arr.U4(1) =  2*V*A^(-1/n)/(xg*ds)^(1+1/n) * (...
        h(2)*((1/n-1)*abs(ud1)^(1/n-2)*ud1+abs(ud1)^(1/n-1)));
    %--------------------------------------------------------------------------
    arr.R(N+1) = -(1+1/n)*2*V*A^(-1/n)/(xg*(xg*ds)^(1+1/n)) * (...
        h(2)*abs(ud1)^(1/n-1)*ud1-h(1)*abs(ud2)^(1/n-1)*ud2)...
        + W*rhog*(h(1)+h(2))*(...
        (1/(2*xg^2*ds))*(h(2)-b(2)-h(1)+b(1))...
        -(1/(2*xg*ds))*(-db(2)+db(1)));
    
    %==========================================================================
    % Interior
    %==========================================================================
    i = 3:N-1;
    %----------------------------------------------------------------------
    % mass conservation
    %----------------------------------------------------------------------
    um1         =  (u(i)  +u(i-1))/2;
    um2         =  (u(i-1)+u(i-2))/2;
    arr.D1(i)   =  -(1/xg).*(um1)./(ds);
    arr.L1(i)   =  -(1/xg).*(-(um2))./(ds);
    %----------------------------------------------------------------------
    % indeed here one equation contains three speeds
    
    arr.D2(i)   =  -(1/xg).*(h(i)*(1/2))./(ds);
    arr.L2(i)   =  -(1/xg).*(h(i)*(1/2)-h(i-1)*(1/2))./(ds);
    arr.LL2(i)  =  -(1/xg).*(-h(i-1)*(1/2))./(ds);
    %----------------------------------------------------------------------
    arr.R(i)    =  (1/xg^2).*(h(i).*(um1)-h(i-1).*(um2))./(ds);
    i = 2:N-1;
    %----------------------------------------------------------------------
    % momentum conservation
    %----------------------------------------------------------------------
    ud1 = u(i+1)-u(i);
    ud2 = u(i)-u(i-1);
    
    arr.D3(i) = 2.*V.*A.^(-1/n)/(xg.*ds).^(1+1/n) .* (-abs(ud2).^(1/n-1).*ud2)...
        - W*rhog.*(1./(2.*xg.*ds)).*(h(i+1)-b(i+1)-h(i)+b(i))...
        + W*rhog.*(h(i)+h(i+1)).*(1./(2.*xg.*ds));
    arr.U3(i) = 2.*V.*A.^(-1/n)./(xg.*ds).^(1+1/n) .* (abs(ud1).^(1/n-1).*ud1)...
        - W*rhog.*(1./(2.*xg.*ds)).*(h(i+1)-b(i+1)-h(i)+b(i))...
        - W*rhog.*(h(i)+h(i+1)).*(1./(2.*xg.*ds));
    %----------------------------------------------------------------------
    arr.D4(i) = 2.*V.*A.^(-1/n)./(xg.*ds).^(1+1/n) .* (...
        -h(i+1).*((1/n-1).*abs(ud1).^(1/n-2).*ud1+abs(ud1).^(1/n-1))...
        -  h(i).*((1/n-1).*abs(ud2).^(1/n-2).*ud2+abs(ud2).^(1/n-1)))...
        -C.*((m-1).*abs(u(i)).^(m-2).*u(i)+abs(u(i)).^(m-1));
    arr.U4(i) = 2.*V.*A.^(-1/n)./(xg.*ds).^(1+1/n) .* (...
        h(i+1).*((1/n-1).*abs(ud1).^(1/n-2).*ud1+abs(ud1).^(1/n-1)));
    arr.L4(i) = 2.*V.*A.^(-1/n)./(xg.*ds).^(1+1/n) .* (...
        h(i).*((1/n-1).*abs(ud2).^(1/n-2).*ud2+abs(ud2).^(1/n-1)));
    %----------------------------------------------------------------------
    arr.R(N+i) = -(1+1/n).*2.*V.*A.^(-1/n)./(xg.*(xg.*ds).^(1+1/n)) .* (...
        h(i+1).*abs(ud1).^(1/n-1).*ud1-h(i).*abs(ud2).^(1/n-1).*ud2)...
        + W*rhog.*(h(i)+h(i+1)).*(...
        (1./(2.*xg.^2.*ds)).*(h(i+1)-b(i+1)-h(i)+b(i))...
        -(1./(2.*xg.*ds))*(-db(i+1)+db(i)));
    %==========================================================================
    % arr.Right boundary
    %==========================================================================
    % mass conservation--------------------------------------------------------
    um1   =  (u(N)  +u(N-1))/2;
    um2   =  (u(N-1)+u(N-2))/2;
    arr.D1(N) =  -(1/xg)*(um1)/(ds);
    arr.L1(N) =  -(1/xg)*(-(um2))/(ds);
    %--------------------------------------------------------------------------
    arr.D2(N)  =  -(1/xg)*(h(N)*(1/2))/(ds);
    arr.L2(N)  =  -(1/xg)*(h(N)*(1/2)-h(N-1)*(1/2))/(ds);
    arr.LL2(N) =  -(1/xg)*(-h(N-1)*(1/2))/(ds);
    %--------------------------------------------------------------------------
    arr.R(N)   =  (1/xg^2)*(h(N)*(um1)-h(N-1)*(um2))/(ds);
    % momentum conservation----------------------------------------------------
    ud2    =  u(N)-u(N-1);
    arr.D4(N)  =  -2*Vb*A^(-1/n)/(xg*ds)^(1/n)*(...
        (1/n-1)*abs(ud2)^(1/n-2)*(ud2) + abs(ud2)^(1/n-1));
    arr.L4(N)  =  -2*Vb*A^(-1/n)/(xg*ds)^(1/n)*(...
        -(1/n-1)*abs(ud2)^(1/n-2)*(ud2) - abs(ud2)^(1/n-1));
    %--------------------------------------------------------------------------
    arr.R(N+N) =   (1/n)*2*Vb*A^(-1/n)/(xg^(1+1/n)*(ds)^(1/n))*...
        abs(u(N)-u(N-1))^(1/n-1)*(u(N)-u(N-1))...
        +Wb*(btrs/2)*(1-rho_ice/rho_water)*rho_water*g*(db(N)+dSL);
    %==========================================================================
    % Auxiliary
    %==========================================================================
    % extrapolation of the flotation criterium
    arr.B(N+N+1) = -2*(rho_water/rho_ice)*(db(N)+dSL);
    arr.B(N)     = +3;
    arr.B(N-1)   = -1;
    
elseif strcmp(model,'roothmod')
    %% Rooth
    
    T1=sol(1);
    T2=sol(2);
    T3=sol(3);
    S1=sol(4);
    S3=sol(5);
    
    k=200;
    
    [V,~,~,~,~,~,Stot] = setvalues(pars);
    
    arr.JR(1,:) = [(T1 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (T1 - T2)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3 - T1 + T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (S1 - S3 - T1 + T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2 - 1,                                                                   -(S1 - S3 - T1 + T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1),                                                             (T1 - T2)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (T1 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (S1 - S3 - T1 + T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 + (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2,                                                                                                                                                            (T1 - T2)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (T1 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 + (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2,                                                                                                                                                            (T1 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (T1 - T2)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2];
    arr.JR(2,:) = [                                        (S1 - S3 - T1 + T3)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) - (T2 - T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (T1 - T2)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) + (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)^2) - (k*exp(-k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2), (S1 - S3 - T1 + T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (S1 - S3 - T1 + T3)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) - 1,                                        (T1 - T2)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) + (T2 - T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (S1 - S3 - T1 + T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)^2) + (k*exp(-k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2),                                                                                                                                            (T1 - T2)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) + (T2 - T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)^2) + (k*exp(-k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2),                                                                                                                                            (k*exp(k*(S1 - S3 - T1 + T3))*(T1 - T2)*(S1 - S3 - T1 + T3))/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)^2) - (T2 - T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (T1 - T2)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) - (k*exp(-k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2)];
    arr.JR(3,:) = [                                                           (T1 - T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (T2 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3 - T1 + T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 + (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2,                                                                     (S1 - S3 - T1 + T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1), (T2 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (T1 - T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3 - T1 + T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (S1 - S3 - T1 + T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (k*exp(k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2 - 1,                                                                                                                                                            (T2 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (T1 - T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (k*exp(k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2,                                                                                                                                                            (T1 - T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (T2 - T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(k*(S1 - S3 - T1 + T3))*(T2 - T3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 + (k*exp(-k*(S1 - S3 - T1 + T3))*(T1 - T3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2];
    arr.JR(4,:) = [                                                                                   (S1 - S3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (S1 + (S1 + S3 - Stot)/V)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (k*exp(k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2,                                                                                                                        0,                                                                                    (S1 + (S1 + S3 - Stot)/V)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 + (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2, (S1 + (S1 + S3 - Stot)/V)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3 - T1 + T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + ((1/V + 1)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 + (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2,         (S1 - S3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (S1 - S3 - T1 + T3)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (S1 + (S1 + S3 - Stot)/V)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) + (S1 - S3 - T1 + T3)/(V*(exp(-k*(S1 - S3 - T1 + T3)) - 1)) - (k*exp(k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2];
    arr.JR(5,:) = [                                                                                   (S1 - S3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) + (S3 + (S1 + S3 - Stot)/V)/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(k*(S1 - S3 - T1 + T3))*(S3 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2,                                                                                                                        0,                                                                                    (k*exp(k*(S1 - S3 - T1 + T3))*(S3 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (S3 + (S1 + S3 - Stot)/V)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2 - (S1 - S3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1),         (k*exp(k*(S1 - S3 - T1 + T3))*(S3 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2 - (S1 - S3 - T1 + T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) - (S3 + (S1 + S3 - Stot)/V)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - (S1 - S3 - T1 + T3)/(V*(exp(k*(S1 - S3 - T1 + T3)) - 1)) - (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2 - (S1 - S3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1), (S1 - S3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) + (S1 - S3 - T1 + T3)/(exp(-k*(S1 - S3 - T1 + T3)) - 1) + (S3 + (S1 + S3 - Stot)/V)/(exp(k*(S1 - S3 - T1 + T3)) - 1) - ((1/V + 1)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1) + (k*exp(-k*(S1 - S3 - T1 + T3))*(S1 - S3)*(S1 - S3 - T1 + T3))/(exp(-k*(S1 - S3 - T1 + T3)) - 1)^2 - (k*exp(k*(S1 - S3 - T1 + T3))*(S3 + (S1 + S3 - Stot)/V)*(S1 - S3 - T1 + T3))/(exp(k*(S1 - S3 - T1 + T3)) - 1)^2];
    
    Jac = arr.JR;
    
elseif strcmp(model,'N')
    %% North
    [NN,drN,A0,n] = setvalues(pars);
    rN   =  -1:drN:1;
    rN = rN(length(rN)/2+1:end-1);
    rNend = 1;
    
    h    =  sol(1:NN);
    hend = 0;
    
    [b,bend] = bedN(rN,pars);
    
    [~,ader] = accNHA(h,pars);
    
    for i=2:NN-1
        arr.DN(i) = ader(i) -((A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(h(i)/2 + h(i-1)/2)^(n + 2))/(2*drN) + (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(h(i)/2 + h(i+1)/2)^(n + 2))/(2*drN) + (A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(n + 2)*(h(i)/2 + h(i-1)/2)^(n + 1)*(b(i) - b(i-1) + h(i) - h(i-1)))/(4*drN) + (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(n + 2)*(h(i)/2 + h(i+1)/2)^(n + 1)*(b(i) - b(i+1) + h(i) - h(i+1)))/(4*drN) + (A0*sign(b(i) - b(i-1) + h(i) - h(i-1))*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 2)*(rN(i) + rN(i-1))*(n - 1)*(h(i)/2 + h(i-1)/2)^(n + 2)*(b(i) - b(i-1) + h(i) - h(i-1)))/(2*drN^2) + (A0*sign(b(i) - b(i+1) + h(i) - h(i+1))*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 2)*(rN(i) + rN(i+1))*(n - 1)*(h(i)/2 + h(i+1)/2)^(n + 2)*(b(i) - b(i+1) + h(i) - h(i+1)))/(2*drN^2))/(drN*rN(i));
        arr.UN(i) = ((A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(h(i)/2 + h(i+1)/2)^(n + 2))/(2*drN) - (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(n + 2)*(h(i)/2 + h(i+1)/2)^(n + 1)*(b(i) - b(i+1) + h(i) - h(i+1)))/(4*drN) + (A0*sign(b(i) - b(i+1) + h(i) - h(i+1))*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 2)*(rN(i) + rN(i+1))*(n - 1)*(h(i)/2 + h(i+1)/2)^(n + 2)*(b(i) - b(i+1) + h(i) - h(i+1)))/(2*drN^2))/(drN*rN(i));
        arr.LN(i) = ((A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(h(i)/2 + h(i-1)/2)^(n + 2))/(2*drN) - (A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(n + 2)*(h(i)/2 + h(i-1)/2)^(n + 1)*(b(i) - b(i-1) + h(i) - h(i-1)))/(4*drN) + (A0*sign(b(i) - b(i-1) + h(i) - h(i-1))*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 2)*(rN(i) + rN(i-1))*(n - 1)*(h(i)/2 + h(i-1)/2)^(n + 2)*(b(i) - b(i-1) + h(i) - h(i-1)))/(2*drN^2))/(drN*rN(i));
    end
    
    for i=1
        arr.DN(i) = ader(i) - (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(h(i)/2 + h(i+1)/2)^(n + 2))/(2*drN^2*rN(i)) - (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(n + 2)*(h(i)/2 + h(i+1)/2)^(n + 1)*(b(i) - b(i+1) + h(i) - h(i+1)))/(4*drN^2*rN(i)) - (A0*sign(b(i) - b(i+1) + h(i) - h(i+1))*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 2)*(rN(i) + rN(i+1))*(n - 1)*(h(i)/2 + h(i+1)/2)^(n + 2)*(b(i) - b(i+1) + h(i) - h(i+1)))/(2*drN^3*rN(i));
        arr.UN(i) = (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(h(i)/2 + h(i+1)/2)^(n + 2))/(2*drN^2*rN(i)) - (A0*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 1)*(rN(i) + rN(i+1))*(n + 2)*(h(i)/2 + h(i+1)/2)^(n + 1)*(b(i) - b(i+1) + h(i) - h(i+1)))/(4*drN^2*rN(i)) + (A0*sign(b(i) - b(i+1) + h(i) - h(i+1))*(abs(b(i) - b(i+1) + h(i) - h(i+1))/drN)^(n - 2)*(rN(i) + rN(i+1))*(n - 1)*(h(i)/2 + h(i+1)/2)^(n + 2)*(b(i) - b(i+1) + h(i) - h(i+1)))/(2*drN^3*rN(i));
        arr.LN(i) = 0;
    end
    
    for i=NN
        arr.DN(i) = ader(i) -((A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(h(i)/2 + h(i-1)/2)^(n + 2))/(2*drN) + (A0*(rNend + rN(i))*(abs(b(i) - bend + h(i))/drN)^(n - 1)*(hend/2 + h(i)/2)^(n + 2))/(2*drN) + (A0*(rNend + rN(i))*(abs(b(i) - bend + h(i))/drN)^(n - 1)*(n + 2)*(hend/2 + h(i)/2)^(n + 1)*(b(i) - bend + h(i)))/(4*drN) + (A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(n + 2)*(h(i)/2 + h(i-1)/2)^(n + 1)*(b(i) - b(i-1) + h(i) - h(i-1)))/(4*drN) + (A0*sign(b(i) - bend + h(i))*(rNend + rN(i))*(abs(b(i) - bend + h(i))/drN)^(n - 2)*(n - 1)*(hend/2 + h(i)/2)^(n + 2)*(b(i) - bend + h(i)))/(2*drN^2) + (A0*sign(b(i) - b(i-1) + h(i) - h(i-1))*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 2)*(rN(i) + rN(i-1))*(n - 1)*(h(i)/2 + h(i-1)/2)^(n + 2)*(b(i) - b(i-1) + h(i) - h(i-1)))/(2*drN^2))/(drN*rN(i));
        arr.UN(i) = 0;
        arr.LN(i) = ((A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(h(i)/2 + h(i-1)/2)^(n + 2))/(2*drN) - (A0*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 1)*(rN(i) + rN(i-1))*(n + 2)*(h(i)/2 + h(i-1)/2)^(n + 1)*(b(i) - b(i-1) + h(i) - h(i-1)))/(4*drN) + (A0*sign(b(i) - b(i-1) + h(i) - h(i-1))*(abs(b(i) - b(i-1) + h(i) - h(i-1))/drN)^(n - 2)*(rN(i) + rN(i-1))*(n - 1)*(h(i)/2 + h(i-1)/2)^(n + 2)*(b(i) - b(i-1) + h(i) - h(i-1)))/(2*drN^2))/(drN*rN(i));
    end
    
    Jac = spdiags([arr.LN([2:NN,1]),arr.DN(1:NN),arr.UN([NN,1:NN-1])],-1:1,NN,NN);%HERE I VHANGES J TO JAC
    Jac = full(Jac);
    
elseif strcmp(model,'SRN')
    %% Full system with full coupling
    [solS,parsS,solR,parsR,solN,parsN] = cut(sol,pars,model);
    
    [NS]=setvalues(parsS);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    [NN,~,~,~,aN,~,~,dadtN,~,DeltaTN] = setvalues(parsN);
    
    TS  = solR(3);
    [parsS(3),der] = Am(TS);
    
    [~]=J(solS,parsS,'S');
    [~]=J(solR,parsR,'roothmod');
    [~]=J(solN,parsN,'N');
    
    % coupling for ocean
    uxg = solS(2*NS);
    hxg = solS(NS);
    %first derivatives wrt hg
    arr.CRS(5,NS) = -sigmaS*uxg;
    %then the derivatives wrt uxg
    arr.CRS(5,2*NS) = -(hxg*sigmaS);
    
    % now for the north
    
    i=1:1:NN;
    [l,dl] = lenN(solN,parsN);
    arr.CRN(4,i)  = - sigmaN * pi* (aN+dadtN*DeltaTN)*2*l*dl(i); 
    
    %coupling for the ice sheet : a new column for the speed equations.
    h    =  sol(1:NS);
    u    =  sol(NS+1:2*NS);
    xg   =  sol(2*NS+1);
    
    [~,dsS,~,~,~,~,~,VS,~,~,n,~,~,~,~,VbS] = setvalues(parsS);
    ud1 = u(2)-u(1);
    ud2 = u(1)+u(1);
    arr.CSR(NS+1) = 2*VS*(-1)/n*parsS(3).^(-1/n-1).*der/(xg*dsS)^(1+1/n) * (...
        h(2)*abs(ud1)^(1/n-1)*ud1-h(1)*abs(ud2)^(1/n-1)*ud2);
    
    i = 2:NS-1;
    ud1 = u(i+1)-u(i);
    ud2 = u(i)-u(i-1);
    
    arr.CSR(NS+i) = 2.*VS.*(-1)/n*parsS(3).^(-1/n-1).*der./(xg.*dsS).^(1+1/n) .* (...
        h(i+1).*abs(ud1).^(1/n-1).*ud1-h(i).*abs(ud2).^(1/n-1).*ud2);
    
    arr.CSR(NS+NS) = -2*VbS*(-1)/n*parsS(3).^(-1/n-1).*der/(xg*dsS)^(1/n)*abs(u(NS)-u(NS-1))^(1/n-1)*(u(NS)-u(NS-1));
    
    %timescaling
    t0mis = t0;
    t0Ro   = t0R;
    
    arr.L1     = arr.L1*t0Ro/t0mis;
    arr.D1     = arr.D1*t0Ro/t0mis;
    arr.U1     = arr.U1*t0Ro/t0mis;
    arr.LL2    = arr.LL2*t0Ro/t0mis;
    arr.L2     = arr.L2*t0Ro/t0mis;
    arr.D2     = arr.D2*t0Ro/t0mis;
    arr.U2     = arr.U2*t0Ro/t0mis;
    arr.R(1:NS) = arr.R(1:NS)*t0Ro/t0mis;
    
elseif strcmp(model,'SRNh')
    %% Full system for hosing
    [solS,parsS,solR,parsR,solN,parsN,newpar] = cut(sol,pars,model);
    
    [NS]=setvalues(parsS);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    [NN,~,~,~,aN,~,~,dadtN,~,DeltaTN] = setvalues(parsN);
    
    TS  = newpar;
    [parsS(3),~] = Am(TS);
    
    [~]=J(solS,parsS,'S');
    [~]=J(solR,parsR,'roothmod');
    [~]=J(solN,parsN,'N');
    
    % coupling for ocean
    uxg = solS(2*NS);
    hxg = solS(NS);
    %first derivatives wrt hg
    arr.CRS(5,NS) = -sigmaS*uxg;
    %then the derivatives wrt uxg
    arr.CRS(5,2*NS) = -(hxg*sigmaS);
    
    % now for the north
    
    i=1:1:NN;
    [l,dl] = lenN(solN,parsN);
    arr.CRN(4,i)  = - sigmaN * pi* (aN+dadtN*DeltaTN)*2*l*dl(i);
    
    %No more coupling in J from the Ocean to South
    
    %timescaling
    t0mis = t0;
    t0Ro   = t0R;
    
    arr.L1     = arr.L1*t0Ro/t0mis;
    arr.D1     = arr.D1*t0Ro/t0mis;
    arr.U1     = arr.U1*t0Ro/t0mis;
    arr.LL2    = arr.LL2*t0Ro/t0mis;
    arr.L2     = arr.L2*t0Ro/t0mis;
    arr.D2     = arr.D2*t0Ro/t0mis;
    arr.U2     = arr.U2*t0Ro/t0mis;
    arr.R(1:NS) = arr.R(1:NS)*t0Ro/t0mis;
    
    
end

end
