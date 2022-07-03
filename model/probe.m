function []=probe(state,statep,pars,dt,model,it)
global prober
% We probe quantities directly in the transient
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

if it==0
else
    
    if strcmp(model,'SRN')
        [stateS,parsS,stateR,~,stateN,parsN] = cut(state,pars,model);
        [statepS,~,~,~,statepN,~] = cut(statep,pars,model);
        
        [~,~,~,~,aN,~,~,dadt,~,DeltaT]=setvalues(parsN);
        [NS] = setvalues(parsS);
        
        [Vol,~]   = VolNmod(stateN,parsN);
        [Volp,~] = VolNmod(statepN,parsN);
        
        xg = stateS(2*NS+1);
        xgp = statepS(2*NS+1);
        
        prober.visc(it) = Am(stateR(3));
        
        prober.outV(it) = -(Vol - Volp)/dt;
        prober.outS(it) = 2*stateS(NS)*stateS(2*NS);
        prober.outSxg(it) = 2*stateS(NS)*(xg - xgp)/dt ;
        prober.outa(it) = pi* (aN+dadt*DeltaT)*lenN(stateN,parsN)^2;
        prober.length(it) = lenN(stateN,parsN);
        prober.accu(it) = aN+dadt*DeltaT;
        prober.VN(it) = Vol;
        prober.VS(it) = VolS(stateS,parsS);
        prober.TN = stateR(1);
        prober.TE = stateR(2);
        prober.TS = stateR(3);
        prober.SN = stateR(4);
        prober.SS = stateR(5);
        
    elseif strcmp(model,'SRNh')
        
        [stateS,parsS,stateR,~,stateN,parsN,newpar] = cut(state,pars,model);
        [statepS,~,~,~,statepN,~] = cut(statep,pars,model);
        
        [~,~,~,~,aN,~,~,dadt,~,DeltaT]=setvalues(parsN);
        [NS] = setvalues(parsS);
        
        [Vol,~]   = VolNmod(stateN,parsN);
        [Volp,~] = VolNmod(statepN,parsN);
        
        xg = stateS(2*NS+1);
        xgp = statepS(2*NS+1);
        
        prober.visc(it) = Am(newpar);
        prober.outV(it) = -(Vol - Volp)/dt;
        prober.outS(it) = 2*stateS(NS)*stateS(2*NS);
        prober.outSxg(it) = 2*stateS(NS)*(xg - xgp)/dt ;
        prober.outa(it) = pi* (aN+dadt*DeltaT)*lenN(stateN,parsN)^2;
        prober.length(it) = lenN(stateN,parsN);
        prober.accu(it) = aN+dadt*DeltaT;
        prober.VN(it) = Vol;
        prober.VS(it) = VolS(stateS,parsS);
        prober.TN = stateR(1);
        prober.TE = stateR(2);
        prober.TS = stateR(3);
        prober.SN = stateR(4);
        prober.SS = stateR(5);
        
    end
    
    
end

end