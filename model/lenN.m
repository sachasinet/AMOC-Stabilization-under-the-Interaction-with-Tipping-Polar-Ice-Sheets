function [l,dl] = lenN (state,pars)
% We compute the adimensional length of Greenland (with theta function)
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

[NN,drN] = setvalues(pars);

%Here is the modification
k=100;
state = 1./(1+exp(-k.*state));

%We compute the len

l = drN/2*state(1);
for i=1:NN
    l = l+state(i)*drN;
end

%We compute the der
dl = zeros(1,NN);
dl(1) = (3*drN*k*exp(-k*state(1)))/(2*(exp(-k*state(1)) + 1)^2);
for i=2:NN
dl(i) = (drN*k*exp(-k*state(i)))/(exp(-k*state(i)) + 1)^2;
end

end