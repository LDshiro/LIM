
function Total=Total_R(coil)
    
Total=0.0;
    
    for i = 1:size(coil,1)
        ss=R_Cu(coil(i,1)*2*pi,coil(i,3)*0.5,20);
        Total=Total+ss;
    end

end
