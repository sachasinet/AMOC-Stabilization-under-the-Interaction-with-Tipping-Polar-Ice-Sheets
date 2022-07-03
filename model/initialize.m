function [initvec,initpars] = initialize(vec,pars,model)
% Here I initilize the system such that the rooth model is in a Northern
% sinking configuration
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

global T0 S0 kdim alphadim betadim topsu tosv T_oceaneq fractoAtl

if strcmp(model,'SRN')
    %% For full Rooth
    [solS,parsS,solR,parsR,solN,parsN] = cut(vec,pars,model);
    
    [NS] = setvalues(parsS);
    [~,~,~,~,aN,~,~,dadtN,~,DeltaTN] = setvalues(parsN);
    [~,~,~,~,~,~,~,sigmaN,sigmaS] = setvalues(parsR);
    
    Rrelax = relax(solR,parsR,'roothmod');
    TSrelax = Rrelax(3);
    T_oceaneq = TSrelax*T0;
    
    parsS(3)=Am(TSrelax); %Not really needed as c=0 at initialization - so only determined by A0S
    Srelax = relax(solS,parsS,'S');
    hxgrelax = Srelax(NS);
    uxgrelax = Srelax(2*NS);
    outS = hxgrelax*sigmaS*uxgrelax;
    parsR(6) = parsR(6)-outS;
    
    Nrelax = relax(solN,parsN,'N');
    outN = sigmaN * pi* (aN+dadtN*DeltaTN)*lenN(Nrelax,parsN)^2;
    parsR(5) = parsR(5)-outN;
    
    initvec = [Srelax',Rrelax',Nrelax']';
    initpars = [parsS,parsR,parsN];
    
    initvec = relax(initvec,initpars,'SRN');
    
    [solS,~,solR,parsR,solN,parsN] = cut(initvec,initpars,model);
    
    TN = solR(1)'*T0;
    TE = solR(2)'*T0;
    TS = solR(3)'*T0;
    SN = solR(4)'*S0;
    SS = solR(5)'*S0;
    
    outN = sigmaN * pi* (aN+dadtN*DeltaTN)*lenN(solN,parsN)^2;
    hxgrelax = solS(NS);
    uxgrelax = solS(2*NS);
    outS = hxgrelax*sigmaS*uxgrelax;
    
    fprintf('\n-- Initialization done --\n')
    fprintf('equilibrium Amoc strength: %2.2e\n',(alphadim*(TS(end)-TN(end))+betadim*(SN(end)-SS(end)))*kdim)
    fprintf('Temperatures: T1 = %2.1f T2 = %2.1f T3 = %2.1f [°C]\n',TN(end),TE(end),TS(end))
    fprintf('Salinities: S1 = %2.1f S3 = %2.1f [psu]\n',SN(end),SS(end))
    fprintf('Outflux: N = %2.3e S = %2.3e [Sv]\n',outN*tosv,outS*tosv/fractoAtl)
    fprintf('Precipitation and melt: N = %2.3e S = %2.3e [psu/s]\n',(outN+parsR(5))*topsu,(outS+parsR(6))*topsu)
    fprintf('Precipitation alone: N = %2.3e S = %2.3e [Sv]\n',parsR(5)*tosv,parsR(6)*tosv)
    fprintf('Viscosity parameter: AS = %2.3e [s-1Pa-3]\n',Am(TS))
    
end
end