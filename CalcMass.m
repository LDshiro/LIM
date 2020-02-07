%CalcMass

%外形、内径から円筒の断面積を計算、→掃引して体積算出→密度から質量変換
function Mass=CalcMass(ra,rb,h)

% ra=50*1e-3;
% rb=0*1e-3;
% h=10e-3;

A=ra^2*pi-rb^2*pi;
V=A*h*1000000;  %m^3→cm^3

Mass=2.70*V*1e-3;   %キログラム

end

