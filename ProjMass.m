%Projectle data gene

function ProjMass=ProjMass(Proj)

ProjMass=0;

 for i = 1:size(Proj,1)
       ProjMass=ProjMass+CalcMass(Proj(i,1)+Proj(i,3)*0.5,Proj(i,1)-Proj(i,3)*0.5,Proj(i,3));  %SelfInd(ä¬îºåa,ÉèÉCÉÑîºåa);
 end

end