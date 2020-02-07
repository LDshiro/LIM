%Cupper wire resistance

function R = R_Cu(length,dia,temp)

%     temp=20; %���x�P���r��
%     length=5.9; %�����@���[�g��
%     dia=0.5*1e-3; %���a ���[�g��

    Cu_roh=(1+0.003862*temp)*1.68e-8;
    %Cu_roh=(1)*1.68e-8;
    R=Cu_roh*length/(dia*dia*pi);

end
