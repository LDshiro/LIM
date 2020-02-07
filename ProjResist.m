%Projectle data gene

function ProjResistMap=ProjResist(Proj)

ProjResistMap=zeros(size(Proj,1),1);

 for i = 1:size(Proj,1)
       ProjResistMap(i)=R_Al(Proj(i,1)*2*pi,Proj(i,3),Proj(i,3),20);  %SelfInd(ä¬îºåa,ÉèÉCÉÑîºåa);
 end                    %R_Al(Proj(i*lpn+j+1,1)^2*pi,ap,ap,20);
                        %R_Al(length,H,W,temp)

end