function Mutual=Mutual_p2c_vec(Coil,Proj)
%Mutual_p2c���x�N�g���ŕԂ���@����ĂԂ��A�������悩�낤�Ǝv���ĂȁB

Mutual=zeros(Proj.NumOfElem,1);

for j=1:Proj.NumOfElem
    
    for i=1:Coil.NumOfTurn
        
        Mutual(j)=Mutual(j)+MutualInd(Coil.CoilMap(i,1),Proj.ElemMap(j,1),Coil.CoilMap(i,2)-Proj.ElemMap(j,2)+Coil.Zc-Proj.Zp);
    
    end
    
end

end
