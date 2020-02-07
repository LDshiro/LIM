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
            %位置マップと抵抗マップ、自己インダクタンスマップは別々に配置しよう。
      end
      
   
   end
end