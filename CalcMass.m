%CalcMass

%O`Aàa©ç~ÌfÊÏðvZA¨|øµÄÌÏZo¨§x©ç¿ÊÏ·
function Mass=CalcMass(ra,rb,h)

% ra=50*1e-3;
% rb=0*1e-3;
% h=10e-3;

A=ra^2*pi-rb^2*pi;
V=A*h*1000000;  %m^3¨cm^3

Mass=2.70*V*1e-3;   %LO

end

