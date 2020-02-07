clear

zc_w=0;             %最下層コイルの最手前座標
rp=20*1e-3;         %プロジェクタイル外径
ac=1*1e-3;          %最下層コイル巻き線太さ
ap=1*1e-3;          %最下層コイル巻き線太さ

rb=1*1e-3;          %バレル厚
rc=rp+rb+ac;        %最下層コイル半径
Num_layC=4;          %Number of layer
Num_lay=1;
lpn=35;             %Number of lay per num

lpnC=17;
z_start=0.0120;

%CoilGene機能に、コイルパラメータ配列を出力する機能を付加しましょう。

CP=CircuitParam;

CP.Rsw=50*1e-3;
CP.Lsw=0.2*1e-6;
CP.RD=20*1e-3;
CP.LD=0.1*1e-6;
CP.C=45*1e-6;
CP.Rci=50*1e-3;
CP.Vc_init=5000;

tspan = [0 0.0004];
options = odeset('RelTol',1e-3,'AbsTol',1e-4);
Coil=RectCoil(rp,rb,ac,ap,Num_layC,lpnC);
Proj=RectProj(rp,ap,Num_lay,lpn);


AA=zeros(Proj.NumOfElem+2+3);
bb=zeros(Proj.NumOfElem+2+3);


%%A行列の定数の入力
for i=3:Proj.NumOfElem+2

    AA(i,i)=Proj.SelfIndMap(i-2);
    bb(i,i)=-Proj.ResistMap(i-2);
    
end

AA(1,1)=CP.Lsw+Coil.SelfInd;
AA(1,2)=CP.Lsw;
AA(2,1)=CP.Lsw;
AA(2,2)=CP.Lsw+CP.LD;
AA(1+2+Proj.NumOfElem,1+2+Proj.NumOfElem)=1;
AA(2+2+Proj.NumOfElem,2+2+Proj.NumOfElem)=1;
AA(3+2+Proj.NumOfElem,3+2+Proj.NumOfElem)=1;


%%これ、定数じゃね？→定数ですわこれ
for i=1:Proj.NumOfElem

    for j=1:Proj.NumOfElem
        
        if i~=j
          AA(2+i,2+j)=MutualInd(Proj.ElemMap(i,1),Proj.ElemMap(j,1),Proj.ElemMap(i,2)-Proj.ElemMap(j,2));
        end
        
    end
    
end


bb(1,1)=-(CP.Rsw+CP.Rci+Coil.Resist);
bb(1,2)=-CP.Rsw;
bb(2,1)=-CP.Rsw;
bb(2,2)=-CP.Rsw-CP.RD;

bb(2+1+Proj.NumOfElem,1)=-1/CP.C;
bb(2+1+Proj.NumOfElem,2)=-1/CP.C;

bb(1,2+1+Proj.NumOfElem)=1;
bb(2,2+1+Proj.NumOfElem)=1;

bb(3+2+Proj.NumOfElem,2+2+Proj.NumOfElem)=1;


%試験的に書いてる
    Proj.Zp=z_start;

y0=zeros(1,Proj.NumOfElem+2+3);
y0(1)=1;

y0(Proj.NumOfElem+2+3)=Proj.Zp;
y0(Proj.NumOfElem+2+1)=CP.Vc_init;

%PCplot
[t,y] = ode45(@(t,y) odeFuncBB_fast(t,y,CP,Coil,Proj,AA,bb), tspan, y0,options);

%plot(t,y(:,Proj.NumOfElem+2+1))
ploter;

vel=y(end,2+Proj.NumOfElem+2);
vel_max=max(y(:,2+Proj.NumOfElem+2));
output_J=0.5*Proj.Mass*vel^2;
input_J=0.5*CP.C*CP.Vc_init^2;
Eff=100*output_J/input_J;




