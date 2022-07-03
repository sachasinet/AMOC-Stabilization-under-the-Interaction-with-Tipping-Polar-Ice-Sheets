function [FEuler]=Euler(F,Fp,theta)
% Applies the theta method on equations
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

%First we modify the rhs F
FEuler = theta*F+(1-theta)*Fp;

%Then we modify the Jacobian entries
global arr
F = fieldnames(arr);
    for iF = 1:length(F)
        aF = F{iF};
        arr.(aF) = arr.(aF) * theta;
    end

end