%Coil data gene

%矩形断面コイルデータ作成
function Coil=CoilGene(rp,rb,ac,ap,Num_lay,lpn)


%zc_w=0;             %最下層コイルの最手前座標
%rp=20*1e-3;         %プロジェクタイル外径
%ac=1*1e-3;          %最下層コイル巻き線太さ
%rb=1*1e-3;          %バレル厚
rc=rp+rb+ac*0.5+ap*0.5;        %最下層コイル半径
% Total_turn=0;

%Num_lay=4;          %Number of layer
%lpn=10;             %Number of lay per num

%コイルに必要なデータ取り出し方
%コイル巻き数をNcとすると、iをNcまでふった時に、すべてのコイルにアクセスできる必要がある。
%つまり、コイルの指定は、単一のアドレスで引き出す必要がある
%必要なデータは、半径データ、ワイヤ径、z座標である。
%z座標は、zc_wからの相対位置が示されていればよい。
%半径は、ワイヤ径を考慮して、矛盾なきよう。
%コイルデータをプロットして、矛盾がないか調べてみよう。

Coil=zeros(Num_lay*lpn,3);

for i=0:Num_lay-1

    for j=0:lpn-1
    
       Coil(i*lpn+j+1,1)=rc+ac*i;
       Coil(i*lpn+j+1,2)=ac*j;      %この値にzc_wを加算すると、ワールド座標になる。
       Coil(i*lpn+j+1,3)=ac;

    end
end

end

