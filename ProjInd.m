%Projectle data gene

function SelfIndProj=ProjInd(Proj)

SelfIndProj=zeros(size(Proj,1),1);

 for i = 1:size(Proj,1)
       SelfIndProj(i)=SelfIndProj(i)+SelfInd(Proj(i,1),Proj(i,3)*0.5);%SelfInd(R,a)
 end

end