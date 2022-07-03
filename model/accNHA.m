function [a,derh,derT] = accNHA(state,parsN,i)
% this is the modified (smooth) accumulation function for the Greenland
% accumulation
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

[~,~,~,~,aN,m,hel,dadt,dhdt,DeltaT] = setvalues(parsN);

k=300;

if nargin < 3
    
    h=state; 
else
    h    =  state(i);
end

a = m.*((h-((aN+dadt.*DeltaT)./m+(hel+dhdt.*DeltaT))) - log(exp((h-((aN+dadt.*DeltaT)./m+(hel+dhdt.*DeltaT))).*k) + 1)./k)+(aN+dadt.*DeltaT);

derh = -m.*(exp(-k.*(hel - h + DeltaT.*dhdt + (aN + DeltaT.*dadt)./m))./(exp(-k.*(hel - h + DeltaT.*dhdt + (aN + DeltaT.*dadt)./m)) + 1) - 1);
derT = dadt - m.*(dhdt + dadt./m - (exp(-k.*(hel - h + DeltaT.*dhdt + (aN + DeltaT.*dadt)./m)).*(dhdt + dadt./m))./(exp(-k.*(hel - h + DeltaT.*dhdt + (aN + DeltaT.*dadt)./m)) + 1));

end