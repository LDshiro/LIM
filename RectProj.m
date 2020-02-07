classdef RectProj
   properties
      NumOfElem
      ElemMap
      ResistMap
      SelfIndMap
      Mass
      Zp=0;
   end
   methods
      function obj = RectProj(rp,ac,Num_lay,lpn)
            obj.ElemMap = ProjGene(rp,ac,Num_lay,lpn);
            obj.SelfIndMap = ProjInd(obj.ElemMap);
            obj.ResistMap=ProjResist(obj.ElemMap);
            obj.Mass=ProjMass(obj.ElemMap);
            obj.NumOfElem=lpn*Num_lay;
            %obj.Zc=zc;
            %�ʒu�}�b�v�ƒ�R�}�b�v�A���ȃC���_�N�^���X�}�b�v�͕ʁX�ɔz�u���悤�B
      end
      
   
   end
end