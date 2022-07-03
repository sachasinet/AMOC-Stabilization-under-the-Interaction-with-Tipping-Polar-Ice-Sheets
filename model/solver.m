function [x,J,P] = solver(F,pars,model)
%======================================================================
% Perform bordered solve on system, return solution.
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
%======================================================================
global arr
if strcmp(model,'S')
% Modified by Sacha Sinet, from Eric Mulder
    %% South
    N = setvalues(pars); % obtain grid size from parameter array
    
    %======================================================================
    % Assemble. I note that this is not the full jacobian as it does not
    % include the flotation equation, nor the derivatives wrt xg
    %======================================================================
    J11 = spdiags([arr.L1([2:N,1]),arr.D1,arr.U1([N,1:N-1])],-1:1,N,N);
    J12 = spdiags([arr.LL2([3:N,1:2]),arr.L2([2:N,1]),...
        arr.D2,arr.U2([N,1:N-1])],-2:1,N,N);
    J21 = spdiags([arr.L3([2:N,1]),arr.D3,arr.U3([N,1:N-1])],-1:1,N,N);
    J22 = spdiags([arr.L4([2:N,1]),arr.D4,arr.U4([N,1:N-1])],-1:1,N,N);
    
    % the sparse function is just a way to save some memory
    J   = sparse([J11,J12;J21,J22]); % so this is F_x
    %======================================================================
    % it seems to me that he is partitioning the Jacobian...
    
    b  = arr.B(1:2*N)'; % there for the flotation
    g  = arr.B(2*N+1);
    r  = arr.R;  % there for the grounding line
    f1 = F(1:2*N);
    f2 = F(2*N+1);
    %======================================================================
    % permutation
    P    = speye(2*N); %This is the sparse identity matrix
    pidx = [(1:N)',(N+1:2*N)']'; % this is what will define the permutation
    P    = P(:,pidx); % this is a way to swp the columns of the identity matrix, to create a permutation matrix
    
    % permute things
    Jp  = P'*J*P;
    f1p = P'*f1;
    rp  = P'*r;
    
    z1 = Jp \ f1p; % c est la qu on solve les systemes, \ cest multiplication par l inverse
    z2 = Jp \ rp;
    
    z1 = P*z1; %z Note that it seems like he uses the change of variables to separate the jacobian into parts to make the system easier to solve.
    z2 = P*z2; %y
    %======================================================================
    x2 = (f2-b*z1)/(g-b*z2); % this is because the system was solving Fxy=Fmu, F_xz=-F, so we need to go back to the original variables!
    x1 = z1 - z2*x2;
    x  = [x1;x2];
    
elseif strcmp(model,'roothmod')
    %% Rooth
    J = arr.JR;
    x =  J\F;
    
elseif strcmp(model,'N')
    %% North
    [NN] = setvalues(pars);
    
    Jh = spdiags([arr.LN([2:NN,1]),arr.DN(1:NN),arr.UN([NN,1:NN-1])],-1:1,NN,NN);%HERE I VHANGES J TO JAC
    x = Jh\F;
    
elseif strcmp(model,'SRN')
    %% SRN
    [rvecS,parsS,rvecR,~,rvecN,parsN] = cut(F,pars,model);
    NS = setvalues(parsS);
    NN = setvalues(parsN);
    
    % We construct the jacobian for the South only
    J11 = spdiags([arr.L1([2:NS,1]),arr.D1,arr.U1([NS,1:NS-1])],-1:1,NS,NS);
    J12 = spdiags([arr.LL2([3:NS,1:2]),arr.L2([2:NS,1]),...
        arr.D2,arr.U2([NS,1:NS-1])],-2:1,NS,NS);
    J21 = spdiags([arr.L3([2:NS,1]),arr.D3,arr.U3([NS,1:NS-1])],-1:1,NS,NS);
    J22 = spdiags([arr.L4([2:NS,1]),arr.D4,arr.U4([NS,1:NS-1])],-1:1,NS,NS);
    Js = sparse([[J11,J12;J21,J22],arr.R;arr.B']);
    
    %We construct the jacobian for the coupled south
    null = zeros(2*NS+1,1);
    Jcs = sparse([null,null,arr.CSR,null,null]);
    
    z    =  solver(rvecS,parsS,'S');
    y    = Js\Jcs;
    w    = solver(rvecN,parsN,'N');
    dstateRguess    = (arr.JR-arr.CRS*y)\(rvecR- arr.CRS * z - arr.CRN * w);
    dstateSguess    = z-y*dstateRguess;
    dstateNguess    = w;
    
    %reconstruction
    x = [dstateSguess',dstateRguess',dstateNguess']';
    
elseif strcmp(model,'SRNh')
    %% hosing
    [rvecS,parsS,rvecR,~,rvecN,parsN] = cut(F,pars,model);
    NS = setvalues(parsS);
    NN = setvalues(parsN);
    
    % We construct the jacobian for the South only
    J11 = spdiags([arr.L1([2:NS,1]),arr.D1,arr.U1([NS,1:NS-1])],-1:1,NS,NS);
    J12 = spdiags([arr.LL2([3:NS,1:2]),arr.L2([2:NS,1]),...
        arr.D2,arr.U2([NS,1:NS-1])],-2:1,NS,NS);
    J21 = spdiags([arr.L3([2:NS,1]),arr.D3,arr.U3([NS,1:NS-1])],-1:1,NS,NS);
    J22 = spdiags([arr.L4([2:NS,1]),arr.D4,arr.U4([NS,1:NS-1])],-1:1,NS,NS);
    Js = sparse([[J11,J12;J21,J22],arr.R;arr.B']);
    
    %We construct the jacobian for the coupled south
    null = zeros(2*NS+1,1);
    Jcs = sparse([null,null,arr.CSR,null,null]);
    
    z    =  solver(rvecS,parsS,'S');
    y    = Js\Jcs;
    w    = solver(rvecN,parsN,'N');
    dstateRguess    = (arr.JR-arr.CRS*y)\(rvecR- arr.CRS * z - arr.CRN * w);
    dstateSguess    = z-y*dstateRguess;
    dstateNguess    = w;
    
    %reconstruction
    x = [dstateSguess',dstateRguess',dstateNguess']';
    
end
end
