function [z,dz] = bed(s,xg,nondim)
% Bedrock definition (By Erik Mulder)
   
    global h0 x0 btype

    if nargin < 3 
        nondim = false;
    end

    if nondim == true
        xg = xg*x0;
        x = s*xg;
    else
        x = s*xg;
    end

    if strcmp(btype, 'type1')
        
        a = 720;
        b = 778.5;
        S = 7e5;
        
        z  = -a + b * (x / S); %z is the b variable of the paper
        dz = b * (s / S);        %dz is db/dx
        
    else
        
        a = 729;
        b = -2184.8;
        c = 1031.72;
        d = -151.72;
        S = 7.5e5;
        z   = -(a + b*(x/S).^2 + c*(x/S).^4 + d*(x/S).^6);
        dz  = -(2*b*xg*(s/S).^2 + 4*c*xg^3*(s/S).^4 + 6*d*xg^5*(s/S).^6);
        
    end

    if nondim == true;
        z  = z / h0;
        dz = dz * x0 / h0;
    end
    
end