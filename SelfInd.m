%Self induntance
%Umoto method
%meter
%�ϐ��ւ̓��͂͂��ׂă��[�g��

function Lu = SelfInd(R,a)

    mu=4*pi*1e-7;
    Li=mu*R/4;
    Le=mu*R*(log(8*R/a)-2);

    Lu=Li+Le;

end

