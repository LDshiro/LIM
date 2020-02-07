function Mutual_dx=Mutual_dx_p2c_vec(Coil,Proj)
%dM/dx���o�͂���ƁA��̓d���͂��Z�o����t�F�[�Y�Ō��ʂ��ė��p�ł������Ȃ̂ŁA�������Ă����āB
%�Ă��A���̂悤�ɂ��܂��傤�B

Mutual_dx=zeros(Proj.NumOfElem,1);

for j=1:Proj.NumOfElem
    
    for i=1:Coil.NumOfTurn
        
        Mutual_dx(j)=Mutual_dx(j)+dMdx(Coil.CoilMap(i,1),Proj.ElemMap(j,1),Coil.CoilMap(i,2)-Proj.ElemMap(j,2)+Coil.Zc-Proj.Zp);
        
    end
    
end

end
