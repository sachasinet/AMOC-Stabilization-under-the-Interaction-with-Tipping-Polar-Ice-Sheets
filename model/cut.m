function [out1,pars1,out2,pars2,out3,pars3,newpar] = cut(vec,pars,model)
% Used to separate pars and state vector into the three systems
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

global svs sps svR spR
    
if strcmp(model,'SRN')
    sizevS = svs;
    sizepS = sps;
    
    sizevR = svR;
    sizepR = spR;
    
    out1 = vec(1:sizevS);
    out2 = vec(sizevS+1:sizevS+sizevR);
    out3 = vec(sizevS+sizevR+1:end);
    
    pars1 = pars(1:sizepS);
    pars2 = pars(sizepS+1:sizepS+sizepR);
    pars3 = pars(sizepS+sizepR+1:end);
    
elseif strcmp(model,'SRNh')
    sizevS = svs;
    sizepS = sps;
    
    sizevR = svR;
    sizepR = spR;
    
    out1 = vec(1:sizevS);
    out2 = vec(sizevS+1:sizevS+sizevR);
    out3 = vec(sizevS+sizevR+1:end);
    
    pars1 = pars(1:sizepS);
    pars2 = pars(sizepS+1:sizepS+sizepR);
    pars3 = pars(sizepS+sizepR+1:end-1); % I change from here, knowing that I have a new parameter
    newpar = pars(end);
    
end

end