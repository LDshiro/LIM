
function Total=Total_Length(coil)
    
Total=0.0;
    
    for i = 1:size(coil,1)
        ss=coil(i,1)*2*pi;
        Total=Total+ss;
    end

end
