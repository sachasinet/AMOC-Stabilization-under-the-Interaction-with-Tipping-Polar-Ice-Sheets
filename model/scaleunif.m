% To adimensionalize the system of coupled equations
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

syms V alpha beta k F1 F3 gamma n T1 T2 T3 S1 S2 S3 tau1 tau2 tau3 lambda Tau1 Tau2 Tau3 Fad1 Fad3 Q V rho_water M y0 ug hg S0 sigma dsN dVolNdt acc L

    T=[T1 T2 T3];
    S=[S1 S2 S3];
    F=[F1,0,F3];
    Fad = [Fad1,0,Fad3];
    tau=[tau1 tau2 tau3];
    Tau=[Tau1 Tau2 Tau3];
    
        for i=[1:3]
        H(i)=lambda.*(tau(i)-T(i));
        end
        F(2) = -1/V.*(F(1)+F(3));
    
    q = k*(alpha*(T(3)-T(1))+beta*(S(1)-S(3)));
        
        out(1) = q*(T(2)-T(1))+H(1);
        out(2) = q/V*(T(3)-T(2))+H(2);
        out(3) = q*(T(1)-T(3))+H(3);
        out(4) = q*(S(2)-S(1))-F(1)+(dVolNdt-acc*pi*L^2)*S0*rho_water/M ;
        out(5) = 1/V*q*(S(3)-S(2))-F(2)+1/V*S0*ug*hg*y0*rho_water/M-(dVolNdt-acc*pi*L^2)*S0*rho_water/M*1/V;
        out(6) = q*(S(1)-S(3))-F(3)-S0*ug*hg*y0*rho_water/M;
        
        %scaling of the variables 
        out=subs(out,[T,S],[T.*lambda/(alpha*k),S.*lambda/(beta*k)]);
        out(1:3)=out(1:3)/(lambda/(alpha*k));
        out(4:6)=out(4:6)./(lambda./(beta*k));
        expand(out(1));
        %time scaling
        out = out/lambda;
        
%       % change of parameters
        out=subs(out,tau,Tau.*lambda./(alpha*k));
        out=subs(out,F,Fad.*lambda^2/(beta*k));
        
        expand(out(4))
        expand(out(6))
        
        % adding the coupling from antarctica
        out=subs(out,rho_water,sigma*M*lambda^2/(S0*beta*k*y0));
        expand(out(5:6));
        