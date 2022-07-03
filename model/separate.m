function [stable, unstable]=separate(cid,curve,pars,model)
% function to separate a continuation curve between stable
% and ustable parts
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

k=size(curve);
stable=zeros(1,k(2));
unstable=zeros(1,k(2));

fprintf('Performing stability analysis ')
for i=1:k(1)
    waiting(i,k(1))
    
    
        pars(cid) = curve(i,end);
        reig = real(eig(J(curve(i,1:end-1),pars,model)));
        
        if sqrt(reig.^2)==-reig
            stable(i,:) = curve(i,:);
        else
            unstable(i,:) = curve(i,:);
        end
        
        
end
    
  stable(find(all(stable==0,2)),:)=[];
  unstable(find(all(unstable==0,2)),:)=[];
  
  fprintf('Stability analysis done \n\n')

end