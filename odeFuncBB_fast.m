function dydt = odeFuncBB_fast(t,y,CP,Coil,Proj,AA,bb)
    
    %dydt = zeros(3,1);

    Proj.Zp=y(2+Proj.NumOfElem+3);
    Ec=y(2+Proj.NumOfElem+1);
    v=y(2+Proj.NumOfElem+2);
    %v=0;
    
    %AA(2,2)=CP.Lsw+CP.LD;
    
    %OdeFuncに書く内容　試験的にこっちに書いてる　ここ、速度が要求されるので、高速な実行を得たい。
    Mpc_vec=Mutual_p2c_vec(Coil,Proj);
    AA(1,3:2+Proj.NumOfElem)=Mpc_vec;
    AA(3:2+Proj.NumOfElem,1)=Mpc_vec;
    %FF=Mutual_p2c_vec(Coil,Proj);

    %試験的に書いてる
    %Coil.Zc=6.9e-3;

    %OdeFuncに書く内容　試験的にこっちに書いてる　ここ、速度が要求されるので、高速な実行を得たい。
    %これ、関数ハンドルうまく使って、値をうまく渡したい。一回だけ実行して変数に書き込む技を身に着けたい。もうグローバル変数でも良い気がするが。

    %vを掛ける
    Mpc_dx_vec=Mutual_dx_p2c_vec(Coil,Proj);
    bb(3:2+Proj.NumOfElem,1)=-v*Mpc_dx_vec;
    bb(1,3:2+Proj.NumOfElem)=-v*Mpc_dx_vec;


    ia=y(1);
    ib=y(2);
    K=ia/Proj.Mass;
    bb(2+2+Proj.NumOfElem,3:2+Proj.NumOfElem)=K*Mpc_dx_vec;

        if Ec>0
            AA(2,2)=CP.Lsw+100;
            %bb(2,2)=-CP.Rsw-CP.RD+100;
            %bb(2+1+Proj.NumOfElem,1)=-1/CP.C;
            %bb(2+1+Proj.NumOfElem,2)=-1/CP.C;
        else
            AA(2,2)=CP.Lsw+CP.LD;
            %bb(2,2)=-CP.Rsw-CP.RD;
            %bb(2+1+Proj.NumOfElem,1)=0;
            %bb(2+1+Proj.NumOfElem,2)=0;
        end

%        AA(2,2)=CP.Lsw+100;
    
        A=inv(AA);

        b1=bb*y;
        
        dydt=A*b1;
    
end
