
function Total_L=Total_Ind(coil)
    
Total_L=0.0;
    
    for i = 1:size(coil,1)
        for j = 1:size(coil,1)
            if i~=j
                ss=MutualInd(coil(i,1),coil(j,1),abs(coil(i,2)-coil(j,2)));
                Total_L=Total_L+ss;
            else
                Total_L=Total_L+SelfInd(coil(i,1),coil(i,3)*0.5);
            end
        end    
    end

end
