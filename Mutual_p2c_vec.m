function Mutual=Mutual_p2c_vec(Coil,Proj)
%Mutual_p2cをベクトルで返すやつ　逐一呼ぶより、効率がよかろうと思ってな。

Mutual=zeros(Proj.NumOfElem,1);

for j=1:Proj.NumOfElem
    
    for i=1:Coil.NumOfTurn
        
        Mutual(j)=Mutual(j)+MutualInd(Coil.CoilMap(i,1),Proj.ElemMap(j,1),Coil.CoilMap(i,2)-Proj.ElemMap(j,2)+Coil.Zc-Proj.Zp);
    
    end
    
end

end
