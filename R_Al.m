%R_Al
function R_Al= R_Al(length,H,W,temp)

    %temp=10+273;　%温度ケルビン
    %length=1000;　%長さ　メートル
    %H  高さ
    %W　幅

    AL_roh=(1+0.0039*temp)*2.82e-8;

    R_Al=AL_roh*length/(H*W);

end