%CalcMass

%�O�`�A���a����~���̒f�ʐς��v�Z�A���|�����đ̐ώZ�o�����x���玿�ʕϊ�
function Mass=CalcMass(ra,rb,h)

% ra=50*1e-3;
% rb=0*1e-3;
% h=10e-3;

A=ra^2*pi-rb^2*pi;
V=A*h*1000000;  %m^3��cm^3

Mass=2.70*V*1e-3;   %�L���O����

end

