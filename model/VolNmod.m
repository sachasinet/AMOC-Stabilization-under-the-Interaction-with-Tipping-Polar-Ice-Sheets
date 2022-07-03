function [out,dstate] = VolNmod (state,pars)
% Computation of the adimensional volume of Greenland (with theta function)
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

[NN,drN] = setvalues(pars);

%Here is the modification - we use a logistic to map negative values to 0.
%Note that it also works at t-1, where h=0 as each
%term is multiplied by h
k=110;
state = state./(1+exp(-k.*state));

rN   =  -1:drN:1;
rN = rN(length(rN)/2+1:end-1);
rNend = 1; 

%We compute the Vol

out = 0;
out = out + 1/2*(state(1)+state(1))*(rN(1)^2-0^2);
for j = 1:NN-1
    out = out + 1/2*(state(j)+state(j+1))*(rN(j+1)^2-rN(j)^2);
end
out = out + 1/2*(state(NN)+0)*(rNend^2-rN(NN)^2);
out = pi*out;

%Now the derivative wrt heights

dstate = zeros(1,NN);

for i = 2:NN-1
dstate(i) = -(pi.*exp(k.*state(i)).*(rN(i-1).^2 - rN(i+1).^2).*(exp(k.*state(i)) + k.*state(i) + 1))./(2.*(exp(k.*state(i)) + 1).^2);
end

i=1;
dstate(i) = (pi*exp(k*state(i))*(rN(i)^2 + rN(i+1)^2)*(exp(k*state(i)) + k*state(i) + 1))/(2*(exp(k*state(i)) + 1)^2);
 
i=NN;
dstate(i) = (pi*exp(k*state(i))*(rNend^2 - rN(i-1)^2)*(exp(k*state(i)) + k*state(i) + 1))/(2*(exp(k*state(i)) + 1)^2);
 

end