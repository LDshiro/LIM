%Cupper wire resistance

function R = R_Cu(length,dia,temp)

%     temp=20; %温度ケルビン
%     length=5.9; %長さ　メートル
%     dia=0.5*1e-3; %半径 メートル

    Cu_roh=(1+0.003862*temp)*1.68e-8;
    %Cu_roh=(1)*1.68e-8;
    R=Cu_roh*length/(dia*dia*pi);

end
