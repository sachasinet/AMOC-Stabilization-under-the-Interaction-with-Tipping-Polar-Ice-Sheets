function []=mainSRN(task)
%==========================================================================
% We run the model containing the WAIS, AMOC and GIS.
% 
%       task: (string) different experiments
% 
%               'STrans': transient for the WAIS alone
%               'RTrans': transient for the AMOC alone
%               'RBif': (pseudo arc length) continuation for the Rooth model
%               'NTrans': transient for the Greenland alone
%               'SRN': transient for the fully coupled model
%               'SRNHos': transient for the hosing experiments including time
%               delays
% 
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
%==========================================================================

global G btype
G=task;
btype = 'type2';

%% Discretization parameters
%used for all (theta method)
theta = 0.7;

%for the South
NS   = 1600;
dsS  =  1/(NS-1); %creates a gri of N points from 0 to 1
sS   =  0:dsS:1;

%for the north

NN=750;
NNN=(NN+1)*2; % to match the half ice sheet symmetrically
drN  =  2/(NNN-1); 
rN   =  -1:drN:1; %this is the total domain, with the two sides
rN = rN(length(rN)/2+1:end-1); %this is the right domain only

%% Dimensional parameters
% used for all 
n         =  3; % Glenn exponent
g         =  9.81; % gravitationnal acceleration
rho_ice   =  910; % Ice density
rho_water =  1000; % Water density
a         =  0.3/(365*24*3600); % present day ice accumulation (from precipitation)

% for the south
AS         =  4.6416e-24;  % Ice viscosity parmeter (redefined in initialization)
betaS      =  0;           % (not used)
CS         =  7.624e6;      % Frictional constant
aS         =  a;
mS         =  1/3;          %1/n
xgS        =  8e5;          % initial grounding line position (not used)
bS         =  bed(sS, xgS, false); % MIS bedrock
btrsS      =  1.0; % buttressing coefficient, 1 stands for no butressing (not used)

rS    = rho_water/rho_ice;

% sea level related (as in Gomez). Later we only consider constant sea
% level.
sl_k   = 4;            % order of fit dS/dxg
sl_Lr  = 20e5;         % location of zero slope
sl_Li  = 6.179253e+05; % approximate location of initial grounding line.

sl_A  = 0;
sl_B  = 0;
sl_C  = 0;

sl = [sl_A, sl_B, sl_C, sl_k, sl_Lr, sl_Li]; % all the parameters are stored in a vector

%for the north (as given by Winkelmann or Greve in General)
LN = 681900; % bedrock radius
mN = 0.005/(365*24*3600); % melting gradient
helN = 1100; % equilibrium line altitude
aN=a; 
AN=1e-16/(365*24*3600); % ice visosity
A0N = 2*AN*(rho_ice*g)^n/(n+2); % rescaled viscosity parameter
ZN = (5*aN*LN^4/(2*AN*(rho_ice*g)^3))^(1/8); %scaling factor

dadtN = 0.015/(365*24*3600); % temprature gradients of a and hel
dhdtN = 100;

DeltaTN = 0; % present day temperature anomaly

%For the Rooth
VR  =  2; % box volume ratio
alphaR =  1.5e-4; % thermal expansion
betaR=8e-4;  %haline expansion
kR=1.5e-6;   %hydraulic constant
F3R=9e-11;   %south precipitation [psu/s]
F1R=3/2*F3R; %north precipitation [psu/s]
tau1R  =  0; %north target temperature 
tau2R  =  30; %equatorial target temperature 
tau3R  =  0; %south target temperature 
lambdaR=  1.29e-9; %newtonian constant for temperature relaxation
MR = 1.08e20; %polar box mass

%% Initial configurations

% South will be defined by initialmax

%for the North
initN = ZN*( 4 * ( (1/2).^ (4/3) - (rN/2).^ (4/3) ) ).^ (3/8);
initN = initN';

%for Rooth
initR = [2.9;28.4;0.3;34.7; 34.1]; % present day as in Lucarini

StotR = 34.7+34.1+VR*35.6;

%% Non-dimensionalization and creation of pars vectors
global h0 x0 u0 t0 t0R T0 S0 svs svn sps spn svR spR kdim alphadim betadim y0 S0R psutosv topsu tosv fractoAtl 

%------------------------------------------------------------------------
% For the South
%------------------------------------------------------------------------
% scales
h0   =  1e3;            % vertical scale in m (south and north)
x0   =  1e5;            % horizontal scale in m (south)
y0 = 2e6/(2*2.795);      % zonal extension south
u0   =  1  / (3600*24);   % typical velocity 1ms^{-1} (south)
t0   =  x0 / u0   ;     % typical time scale (south)

%rescaled pars
aS    =  aS*x0/h0/u0;
rhogS =  rho_ice*g;
VS    =  (h0 / x0) * (1/t0)^(1/n) / CS;
WS    =  (h0^2 / x0) / CS;
VbS   =  (1/t0)^(1/n)/h0;  % at the right bdy
WbS   =  1;
M0S   =  1; %(h0 / t0); not used
CS    =  (u0)^mS;

%pars for the function Am
global T_oceaneq cAm A0S
T_oceaneq = 0.2975; %initial temperature of the southern box
cAm = 0; %coupling constant gamma_S
A0S = 2e-25 ; %initial viscosity parameter

% on cree les vecteurs state et param
parsS = [NS,dsS,AS,betaS,aS,rhogS,CS,VS,WS,mS,n,g,rho_ice,rho_water,sl,VbS, ...
    WbS,M0S,h0,x0,u0,t0,btrsS];

% initS  = [hS;u;xg], this is the used initial state to be relaxed
initS=readmatrix('initialmax.dat');
parsS(3)=initS(2*NS+2); %For the moment at 2.46e-26
initS = initS(1:2*NS+1);

[svs,~] = size(initS);
[~,sps] = size(parsS);

%------------------------------------------------------------------------
% For rooth
%------------------------------------------------------------------------
t0R = 1/lambdaR; %timescale for Rooth and Greenland
kdim = kR*MR/rho_water;
alphadim = alphaR;
betadim = betaR;
S0R = 35; %average salinity of the Atlantic
psutosv = MR/rho_water/S0R*1e-6; %conversion constant
topsu = S0/t0R; %conversion constant
tosv = S0/t0R*MR/rho_water/S0R*1e-6; %conversion constant

T0 = lambdaR/(alphaR*kR); %scales temperature
S0 = lambdaR/(betaR*kR); % scaled salinity

tau1R = tau1R/T0;
tau2R = tau2R/T0;
tau3R = tau3R/T0;
F1R = F1R/S0*t0R;
F3R = F3R/S0*t0R;
fractoAtl = 0.27;

StotR = StotR/S0;

sigmaN = rho_ice/MR * LN^2*h0 * S0R/S0;

sigmaS = 2*fractoAtl*rho_ice/MR * y0*h0*u0 * S0R/S0 * t0R;

parsR = [VR,tau1R,tau2R,tau3R,F1R,F3R,StotR,sigmaN,sigmaS];

initR(1:3) = initR(1:3)/T0;
initR(4:5) = initR(4:5)/S0;

[svR,~] = size(initR);
[~,spR] = size(parsR);

%------------------------------------------------------------------------
%For North
%------------------------------------------------------------------------
%note: directly scaled on the time of Rooth model
A0N = A0N/(LN^(n+1))*h0^(2*n+1)*t0R;
aN = aN/h0*t0R;
mN = mN*t0R;
helN = helN/h0;
dadtN = dadtN/h0*t0R*T0;
dhdtN = dhdtN/h0*T0;
DeltaTN = DeltaTN/T0;

parsN = [NN,drN,A0N,n,aN,mN,helN,dadtN,dhdtN,DeltaTN];

initN = initN/h0; %initial state to be relaxed

[svn,~] = size(initN);
[~,spn] = size(parsN);

%% Allocating Jacobian vectors
%==========================================================================
% Global struct that will contain the Jacobian: for the south
%     #1(i) stands for the derivative wrt h(i) in mass
%     #2(i) stands for the derivative wrt u(i) in mass
%     #R(i) stands for the derivative wrt x_g  in mass
%     #3(i) stands for the derivative wrt h(i) in mom
%     #4(i) stands for the derivative wrt u(i) in mom
%     #R(i) stands for the derivative wrt x_g  in mom
%     #B(i) stands for the derivatives in the flotation
%
%==========================================================================
global arr
arr = [];
arr.D1  = zeros(NS,1);  arr.D2  = zeros(NS,1);
arr.D3  = zeros(NS,1);  arr.D4  = zeros(NS,1);
arr.U1  = zeros(NS,1);  arr.U2  = zeros(NS,1);
arr.U3  = zeros(NS,1);  arr.U4  = zeros(NS,1);
arr.L1  = zeros(NS,1);  arr.L2  = zeros(NS,1);
arr.LL2 = zeros(NS,1);  arr.L3  = zeros(NS,1);
arr.L4  = zeros(NS,1);
arr.R   = zeros(2*NS,1);
arr.B   = zeros(2*NS+1,1);

% And now we add the components for Rooth.
arr.JR   = zeros(5,5);

% And now we add the components for the north.
arr.DN   = zeros(NN,1);
arr.LN   = zeros(NN,1);
arr.UN   = zeros(NN,1);
arr.LLN   = zeros(NN,1);
arr.UUN   = zeros(NN,1);
arr.Glob = zeros(NN,NN);

% coupling from rooth to ice caps
arr.CSR = zeros(2*NS+1,1);
arr.CRS  = zeros(5,2*NS+1);
arr.CRN  = zeros(5,NN);

%% Intialize Global
%Here we define the initial state presented in the Suplementary Material (via the function initialize)
if ~strcmp(task,'RTrans') && ~strcmp(task,'RBif')
pars = [parsS,parsR,parsN];
    init = [initS',initR',initN']';
    [init,pars] = initialize(init,pars,'SRN');
    [initS,parsS,initR,parsR,initN,parsN] = cut(init,pars,'SRN');
end

%% experiments

if strcmp(G,'STrans')
    %% transient South
    %======================================================================
    % Here we perform the transient for the MIS alone
    %======================================================================
    %------------------------------------------
    %First we set the discretization parameters
    %------------------------------------------
    dt = 2; %the adimensional time  
    K  = 150; %number of timesteps
    
    %------------------------------------------
    %Now we set the parameters for the MIS
    %------------------------------------------
    
    prevS = relax(initS,parsS,'S');
    xgplot = zeros(K,1);
    hg = zeros(K,1);
    ug = zeros(K,1);
    vol = zeros(K,1);
    
    for i=1:K
        if i>1
            parsS(3)=2.900e-25; %example
        end
        
        [transS]=transient(prevS,parsS,dt,theta,'S',i);
        prevS = transS;
       
        
        vol(i)=VolS(transS,parsS);
        
        % Here we extract the variables to plot
        h = transS(1:NS);
        hg(i) = transS(NS)*h0;
        ug(i) = transS(2*NS)*u0;
        xgplot(i) = transS(2*NS+1)*x0;
        xg = transS(2*NS+1);
        b=bed(sS, xg , true); %true because s and xg are not dimensional there
        bplot=bed(sS, 16 , true);
        SL=zeros(1,NS);
        
        % and we dimensionalize
        h = h*h0;
        xg=xg*x0;
        b=b*h0;
        bplot=bplot*h0;
        xgplot(i)=xg;
        
        % And we plot
        figure(1)
        nexttile(1)
        plot(sS*xg,h'-b);
        hold on
        plot(sS*16*x0,-bplot)
        plot(sS*16*x0,SL)
        legend('ice','ground','sea')
        hold off
        xlabel('distance [m]')
        ylabel('height [m]')
        ylim([-3000 5000])
        
    end
    
elseif strcmp(G,'RTrans')
    %% transient rooth
    initR = relax(initR,parsR,'roothmod');
    
    dtyrs = 10;
    K=300;
    dt = dtyrs*(365*24*3600)/t0R;
    
    trans = zeros(K,5);
    trans(1,:) = initR;
    
    
    %we build a 1000 yrs temperature forcing
    t00 = 1000;
    t00 = t00/t0R*(365*24*3600);
    tau1ttot = 6.5/T0;
    tau3ttot = 0.3*tau1ttot;
    
    Kstop = t00/dt+1; 
    DF1 = tau1ttot/(Kstop-1);
    DF3 = tau3ttot/(Kstop-1);
    
    previous = initR;
    for i=1:K
        if i>1 && i<=Kstop
            parsR([2,4])=parsR([2,4])+[DF1,DF3];
        end
        
        state = transient(previous,parsR,dt,theta,'roothmod');
        previous=state;
        trans(i,:)=state';
        
    end
    figure(2)
    plot([0:K-1]*dtyrs,((trans(:,3)-trans(:,1))+(trans(:,4)-trans(:,5)))./0.1140)
    title('AMOC strength ')
    xlabel('t [yrs]','Interpreter','latex')
    ylabel('q(t)/q(0)','Interpreter','latex')
    
elseif strcmp(G,'RBif')
    %% bifurcation rooth
    dr = 0.01; %continuation parameter
    dp = 0.01;
    kc = 800;
    parsc = [dr,dp];
    initR = relax(initR,parsR,'roothmod');
    
    curve = continuation(5,initR,parsR,parsc,kc,'roothmod');
    
    [stable, unstable]=separate(5,curve,parsR,'roothmod');
    
    figure(2)
    plot(stable(:,end).*tosv,kR.*((stable(:,3)-stable(:,1)).*T0*alphaR+(stable(:,4)-stable(:,5)).*S0*betaR).*MR./rho_water/1e6,'Linestyle','none','Marker','.','Color','b')
    hold on
    plot(unstable(:,end).*tosv,kR.*((unstable(:,3)-unstable(:,1)).*T0*alphaR+(unstable(:,4)-unstable(:,5)).*S0*betaR).*MR./rho_water/1e6,'Linestyle','none','Marker','.','Color','r')
    hold off
    legend('stable','unstable')
    title('bifurcation diagram')
    ylabel('AMOC strength $q$ [Sv]','Interpreter','latex')
    xlabel('FN [Sv]')
    
elseif strcmp(G,'NTrans')
    %% transient N
    
    figure(1)
    plot(rN.*LN,initN*h0)
    hold on
    
    
    h0N = initN*1.05;
    
    initN = relax(h0N,parsN,'N');
        
    K=200;
    dtyrs = 10;
    dt = dtyrs/t0R*(365*24*3600);
    
    Len = zeros(K,1);
    Len(1) = 1;
    
    prec = zeros(K,1);
    prec(1) = aN;
    
    Vol = zeros(K,1);
    Vol(1) = VolN(initN,parsN);
    
    Volm = zeros(K,1);
    Volm(1) = VolNmod(initN,parsN);
    
    Lenm = zeros(K,1);
    Lenm(1) = 1;
    
    for i = 1:K
        
        if i>1
          parsN(end) = 10/T0; % it is a one shot forcing  
        end
        
        state = transient(initN,parsN,dt,theta,'N');
        
        %In what follows, I compare the Vol and length with or without the
        %use of the theta function. First with the theta function
        Volm(i) = VolNmod(state,parsN);
        [Lenm(i),~] = lenN(state,parsN);
        
        % Now without the theta function
        state = max(state,0);
        initN=state;
        statelen = length(state(state>0));
        Len(i) = drN*statelen+drN/2;
        if statelen==0
            Len(i) = 0;
        end       
        Vol(i) = VolN(state,parsN);
        
        prec(i) = aN + dadtN*parsN(end);
        
        hold on
        if i==1 || mod(i,50)==0
            figure(2)
            plot([rN,1].*LN,[(state*h0)',0])
            xlim([0 750000])
            xlabel('x [km]')
            ylabel('h [m]')
        end
        
        
    end
    hold off
    
    % I construc the Vol and outflux vectors 
    dVol = zeros(K-1,1);
    dVolm = zeros(K-1,1);
    j = 1:K-1;
    dVol(j) = (Vol(j+1)-Vol(j))/(dt);
    dVolm(j) = (Volm(j+1)-Volm(j))/(dt);
    
    figure(3)
    nexttile(1)
    plot([1:K-1].*dtyrs,-dVol.*h0*LN^2/t0R*1e-6)
    hold on
    plot([1:K-1].*dtyrs,-dVolm.*h0*LN^2/t0R*1e-6)
    legend('Original','Smoothened')
    title('Comparison of methods Outfluxes')
    hold off
    
    nexttile(2)
    plot([0:K-1].*dtyrs,Vol*h0*LN^2)
    hold on
    plot([0:K-1].*dtyrs,Volm*h0*LN^2)
    hold off
    legend('Original','Smoothened')
    title('Comparison of methods for Volumes')
    
    nexttile(3)
    plot([0:K-1].*dtyrs,Len.*LN)
    hold on
    plot([0:K-1].*dtyrs,Lenm.*LN)
    hold off
    legend('Original','Smoothened')
    title('Comparison of methods for Lengths')
    
elseif strcmp(G,'SRN')
    %% transient full model with amplification
    %Full coupled model. Here we directly work with amplification coefficients
    
    dtyrs = 10;
    K = 300;
    dt = dtyrs*(365*24*3600)/t0R;
    
    tau2ttot = 5.1/T0; % We set the equatorial temperature change
    %Amplification coefficients
    ampN = 1.3; %Northern box
    ampS = 1;    % Southern box
    
    %coupling constant
    cAm = 0;
    
    global prober
    prober = [];
    prober.outV = zeros(1,K);
    prober.outa = zeros(1,K);
    prober.outS = zeros(1,K);
    prober.length = zeros(1,K);
    prober.VN = zeros(1,K);
    prober.VS = zeros(1,K);
    prober.visc = zeros(1,K);
    prober.outSxg = zeros(1,K);
    
    transS = zeros(K,2*NS+1);
    transS(1,:) = initS;
    
    transR = zeros(K,5);
    transR(1,:) = initR;
    
    transN = zeros(K,NN);
    transN(1,:) = initN;
    
    %We construct the 100 yrs forcing
    t00 = 100;
    t00 = t00*(365*24*3600)/t0R; % we adimentionalize
    
    tau1ttot = ampN*tau2ttot;
    tau3ttot = ampS*tau2ttot;
    
    fprintf('Global average T change: %2.2e\n',computeglob(tau2ttot*T0))
    fprintf('T change in northern Atlantic: %2.2e\n',tau1ttot*T0)
    fprintf('T change in Eq: %2.2e\n',tau2ttot*T0)
    fprintf('T change in southern Atlantic: %2.2e\n',tau3ttot*T0)
    fprintf('T change in North pole: %2.2e\n',tau1ttot/1.25*2.25*T0)
    
    Kstop = t00/dt+1; 
    DF1 = tau1ttot/(Kstop-1);
    DF2 = tau2ttot/(Kstop-1);
    DF3 = tau3ttot/(Kstop-1);
    
    for i=1:K
        if i~=1 && i<=Kstop
            parsR([2,3,4])=parsR([2,3,4])+[DF1,DF2,DF3];
            parsN(10) = parsN(10)+2/ampN*DF1;
            pars = [parsS,parsR,parsN];
        end
        
        state = transient(init,pars,dt,theta,'SRN',i);
        [stateS,~,stateR,~,stateN] = cut(state,pars,'SRN');
        stateN = max(stateN,0);
        init=[stateS',stateR',stateN']';
        [transS(i,:),~,transR(i,:),~,transN(i,:)] = cut(init,pars,'SRN');
        
        DeltaTN(i) = parsN(10);
        
        
    end
    
    %rescaling
    hS = transS(:,1:NS)'*h0;
    u = transS(:,NS+1:2*NS)'*u0;
    xg = transS(:,2*NS+1)'*x0;
    
    hN = transN(:,1:NN)'*h0;
    
    TN = transR(:,1)'*T0;
    TE = transR(:,2)'*T0;
    TS = transR(:,3)'*T0;
    SN = transR(:,4)'*S0;
    SS = transR(:,5)'*S0;
    
    figure(1)
    nexttile(1)
    plot([0:K-1]*dtyrs,kR*(alphaR*(TS-TN)+betaR*(SN-SS))*MR/rho_water./1e6)
    title('AMOC strength')
    xlabel('time [yrs]')
    ylabel('q [Sv]')
    nexttile(2)
    plot([0:K-1]*dtyrs,xg)
    title('grounding line position S')
    xlabel('time [yrs]')
    ylabel('x_g [m]')
    nexttile(3)
    plot([0:K-1]*dtyrs,prober.VN.*h0*LN^2)
    hold on
    plot([0:K-1]*dtyrs,prober.VS.*h0*x0*y0)
    hold off
    title('Ice sheet volumes')
    legend('North','South')
    xlabel('time [yrs]')
    ylabel('V [m^3]')
    nexttile(4)
    plot([0:K-1]*dtyrs,TS)
    hold on
    plot([0:K-1]*dtyrs,TN)
    hold off
    title('Temperatures')
    legend('T south', 'T North')
    xlabel('time [yrs]')
    ylabel('T [°C]')
    nexttile(5)
    plot([0:K-1]*dtyrs,SS)
    hold on
    plot([0:K-1]*dtyrs,SN)
    hold off
    title('Salinities')
    xlabel('time [yrs]')
    ylabel('S [psu]')
    
    nexttile(7)
    plot([0:K-1]*dtyrs,prober.visc)
    title('Viscosity')
    xlabel('time [yrs]')
    ylabel('Viscosity')
    
    fprintf('\n-- trasient done --\n')
    fprintf('equilibrium Amoc strength: %2.2e\n',kdim*(alphadim*(TS(end)-TN(end))+betadim*(SN(end)-SS(end))))
    fprintf('equilibrium temperatures: T1 = %2.1f T2 = %2.1f T3 = %2.1f \n',TN(end),TE(end),TS(end))
    fprintf('equilibrium salinities: S1 = %2.1f S3 = %2.1f \n',SN(end),SS(end))
    
    toplotN = prober.VN./prober.VN(1);
    toplotS = prober.VS./prober.VS(1); 
    
    figure(2)
    plot([0:K-1]*dtyrs,kR*(alphaR*(TS-TN)+betaR*(SN-SS))*MR/rho_water./(kR*(alphaR*(TS(1)-TN(1))+betaR*(SN(1)-SS(1)))*MR/rho_water))
    hold on
    plot([0:K-1]*dtyrs,toplotN)
    plot([0:K-1]*dtyrs,toplotS)
    hold off
    ylim([-4 2])
    xlabel('time [yr]','Interpreter','latex')
    legend('$q$','$V_{\mbox{N}}$','$V_{\mbox{S}}$','Interpreter','latex')
    title('Full transient')
    
elseif strcmp(G,'SRNHos')
    %% Same as SRN but for hosing experiment with time delays
    
    cdelexp = [0,0.175;0,0.1;1100,1;1700,1]; % we set different couplings
    [nexp,~] = size(cdelexp);
    subt = ["(a) No AMOC tipping","(b) AMOC tipping","(c) No AMOC tipping","(d) AMOC tipping"];
    
    parssave = pars;
    initsave=init;
    
    counter=0;
    figure(2)
    t = tiledlayout(2,2);
    
    for exp = 1:nexp
        counter = counter + 1;
        cit = cdelexp(exp,2);
        delayit = cdelexp(exp,1);
        
        init = initsave;
        pars = parssave;
        [initS,parsS,initR,parsR,initN,parsN] = cut(init,pars,'SRN');
        
        newpar = initR(3);
        pars(end+1) = newpar;
        
        dtyrs = 5;
        K = 600;
        dt = dtyrs*(365*24*3600)/t0R;
        
        global prober
        prober = [];
        prober.outV = zeros(1,K);
        prober.out2 = zeros(1,K);
        prober.outa = zeros(1,K);
        prober.outS = zeros(1,K);
        prober.length = zeros(1,K);
        prober.VN = zeros(1,K);
        prober.VS = zeros(1,K);
        prober.visc = zeros(1,K);
        prober.outSxg = zeros(1,K);
        
        transS = zeros(K,2*NS+1);
        transS(1,:) = initS;
        
        transR = zeros(K,5);
        transR(1,:) = initR;
        
        transN = zeros(K,NN);
        transN(1,:) = initN;
        
        
        DelayN = 1500;
        DelayN = DelayN*(365*24*3600)/t0R;
        KbeginN = DelayN/dt+1;
        
        DelayS = delayit;
        DelayS = DelayS*(365*24*3600)/t0R;
        KbeginS = DelayS/dt+1;
        
        tforc = 100;
        tforc = tforc*(365*24*3600)/t0R;
        
        tau1ttot = 23/T0;
        KstopN = tforc/dt+KbeginN;
        DF1 = tau1ttot/(KstopN-KbeginN);
        
        TStot = 7/T0;
        KstopS = tforc/dt+KbeginS;
        DFn = TStot/(KstopS-KbeginS);
        
        cAm = cit;
        
        for i=1:K
            if i~=1 & i<=KstopN & i>KbeginN
                parsN(10) = parsN(10)+DF1;
            end
            
            if i~=1 & i<=KstopS & i>KbeginS
                newpar = newpar + DFn;
            end
            
            pars = [parsS,parsR,parsN,newpar];
            
            state = transient(init,pars,dt,theta,'SRNh',i);
            [stateS,~,stateR,~,stateN] = cut(state,pars,'SRNh');
            stateN = max(stateN,0);
            init=[stateS',stateR',stateN']';
            [transS(i,:),~,transR(i,:),~,transN(i,:)] = cut(init,pars,'SRNh');
            
            DeltaTN(i) = parsN(10);
            
            
        end
        
        TN = transR(:,1)'*T0;
        TE = transR(:,2)'*T0;
        TS = transR(:,3)'*T0;
        SN = transR(:,4)'*S0;
        SS = transR(:,5)'*S0;
        
        q = kdim*(alphadim*(TS(end)-TN(end))+betadim*(SN(end)-SS(end)));
        if q>0
            disp('notip')
        else
            disp('tip')
        end
        
        nexttile
        plot([0:K-1]*dtyrs,prober.outV.*h0*LN^2/t0R.*1e-6*rho_ice/rho_water+prober.outa.*h0*LN^2/t0R.*1e-6*rho_ice/rho_water,'LineWidth',1)
        hold on
        plot([0:K-1]*dtyrs,fractoAtl*prober.outS.*h0*y0*u0.*1e-6*rho_ice/rho_water-fractoAtl*prober.outSxg.*h0.*y0*x0/t0R.*1e-6*rho_ice/rho_water,'LineWidth',1)
        hold off
        subtitle(subt(counter),'Interpreter','latex')
        legend('GIS','WAIS','Interpreter','latex')
        
        if counter == 1
            set(gca, 'xticklabel', [])
        elseif counter ==2
            set(gca, 'xticklabel', [])
            set(gca, 'yticklabel', [])
        elseif counter ==4
            set(gca, 'yticklabel', [])
        end
           
    end
    
    % Add shared title and axis labels
    xlabel(t,'Time [years]','Interpreter','latex')
    ylabel(t,'Hosing [Sv]','Interpreter','latex')
%     % Move plots closer together
    t.TileSpacing = 'compact';

end




