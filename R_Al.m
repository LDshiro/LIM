%R_Al
function R_Al= R_Al(length,H,W,temp)

    %temp=10+273;�@%���x�P���r��
    %length=1000;�@%�����@���[�g��
    %H  ����
    %W�@��

    AL_roh=(1+0.0039*temp)*2.82e-8;

    R_Al=AL_roh*length/(H*W);

end