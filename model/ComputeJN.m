syms NN drN rNi rNend hend rNip rNim dsN A0 n aN hip hipp him himm hi h2 h1 hNNm hNN himmm himmmm himmmmm bi bim bip bend
% Used to compute the jacobian of the Greenland model  
    Hi = bi+hi;
    Him = bim+him;
    Hip = bip+hip;
    Hend = bend;
    %% For the Centre
    Di = A0*(rNi+rNip)/2*((hi+hip)/2)^(n+2)*(abs(Hip-Hi)/(drN))^(n-1);
    Dim = A0*(rNim+rNi)/2*((him+hi)/2)^(n+2)*(abs(Hi-Him)/(drN))^(n-1);
    
    Fluxi = -Di*(Hip-Hi)/(drN);
    Fluxim = -Dim*(Hi-Him)/(drN);
    
    dFluxi = 1/(rNi*drN)*(Fluxi-Fluxim);
    
    outi = -dFluxi+aN;
    
    fprintf('\n')
    fprintf('for i=2:NN-1 \n')
    fprintf('arr.DNm(i) = %s;\n',diff(outi,hi))
    fprintf('arr.UNm(i) = %s;\n',diff(outi,hip))
    fprintf('arr.LNm(i) = %s;\n',diff(outi,him))
    fprintf('end \n\n')
    
    %% For the extreme - Left
    Di = A0*(rNi+rNip)/2*((hi+hip)/2)^(n+2)*(abs(Hip-Hi)/(drN))^(n-1);
    
    Fluxi = -Di*(Hip-Hi)/(drN);
    Fluxim = 0;
    
    dFluxi = 1/(rNi*drN)*(Fluxi-Fluxim);
    
    outi = -dFluxi+aN;
    
    fprintf('for i=1 \n')
    fprintf('arr.DNm(i) = %s;\n',diff(outi,hi))
    fprintf('arr.UNm(i) = %s;\n',diff(outi,hip))
    fprintf('arr.LNm(i) = %s;\n',diff(outi,him))
    fprintf('end \n\n')
    
     %% For the extreme - right
    Di = A0*(rNi+rNend)/2*((hi+hend)/2)^(n+2)*(abs(Hend-Hi)/(drN))^(n-1);
    Dim = A0*(rNim+rNi)/2*((him+hi)/2)^(n+2)*(abs(Hi-Him)/(drN))^(n-1);
    
    Fluxi = -Di*(Hend-Hi)/(drN);
    Fluxim = -Dim*(Hi-Him)/(drN);
    
    dFluxi = 1/(rNi*drN)*(Fluxi-Fluxim);
    
    outi = -dFluxi+aN;
    
    fprintf('for i=NN \n')
    fprintf('arr.DNm(i) = %s;\n',diff(outi,hi))
    fprintf('arr.UNm(i) = %s;\n',diff(outi,hip))
    fprintf('arr.LNm(i) = %s;\n',diff(outi,him))
    fprintf('end \n\n')