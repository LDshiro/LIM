%Projectle data gene

function Proj=ProjGene(rp,ap,Num_lay,lpn)
%Proj(最外半径、分解能、層数、層あたり巻き数)
%矩形断面プロジェクタイルデータ作成

%rp=20*1e-3;         %スリーブ最外径
%ap=1*1e-3;          %スリーブ分解能
%Num_lay=3;          %Number of layer
%lpn=10;             %Number of lay per num

%コイルに必要なデータ取り出し方と同じ
%コイル巻き数をNcとすると、iをNcまでふった時に、すべてのコイルにアクセスできる必要がある。
%つまり、コイルの指定は、単一のアドレスで引き出す必要がある
%必要なデータは、半径データ、ワイヤ径、z座標である。
%z座標は、zc_wからの相対位置が示されていればよい。
%半径は、ワイヤ径を考慮して、矛盾なきよう。
%コイルデータをプロットして、矛盾がないか調べてみよう。

Nc=Num_lay*lpn;

Proj=zeros(Nc,3);

for i=0:Num_lay-1

    for j=0:lpn-1
    
       Proj(i*lpn+j+1,1)=rp-ap*i;   %最外径を基準に、内側に積層する感じ
       Proj(i*lpn+j+1,2)=ap*j;      %この値にzp_wを加算すると、ワールド座標になる。
       Proj(i*lpn+j+1,3)=ap;        %ワイヤ太さ（直径
%        Proj(i*lpn+j+1,4)=R_Al(Proj(i*lpn+j+1,1)^2*pi,ap,ap,20);         %R_Al(length,H,W,temp)
%        Proj(i*lpn+j+1,5)=SelfInd(Proj(i*lpn+j+1,1),Proj(i*lpn+j+1,3)*0.5);  %SelfInd(環半径,ワイヤ半径);
    end
end


end