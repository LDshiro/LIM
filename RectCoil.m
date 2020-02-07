classdef RectCoil
   properties
      NumOfTurn
      WireLength
      SelfInd
      Resist
      CoilMap
      Zc=0;
   end
   methods
%       function CoilMap = CoilGene(rp,ac,Num_lay,lpn)
%          CoilMap = CoilGene(rp,ac,Num_lay,lpn);
%       end
      function obj = RectCoil(rp,rb,ac,ap,Num_lay,lpn)
            obj.CoilMap = CoilGene(rp,rb,ac,ap,Num_lay,lpn);
            obj.Resist=Total_R(obj.CoilMap);
            obj.SelfInd=Total_Ind(obj.CoilMap);
            obj.WireLength=Total_Length(obj.CoilMap);
            obj.NumOfTurn=lpn*Num_lay;
            %obj.Zc=zc;
      end
      
   
   end
end