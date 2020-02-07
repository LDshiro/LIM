%gradienr coil mapper

function Field = BioSava_gradient2(A,zc,zp,FieldSize,I)

Field=zeros(FieldSize,FieldSize);
Field_vec=zeros(FieldSize*FieldSize,1);

cnt=size(A);

%改良の方向としては、Aを連続のベクトルに変換する方向で取り組もう。
%あ?でも、Aは前後の差分を取るベクトルでもあるのでな?
%えっと、

% x=A(j,1,i)*1e-3;
% y=A(j,2,i)*1e-3;
% xd=(A(j,1,i)-A(j+1,1,i))*1e-3;
% yd=(A(j,2,i)-A(j+1,2,i))*1e-3;

%あるiに対して、この情報を持った4列のベクトルを用意してあげれば、面倒な計算を一気に畳込める。

%評価領域もベクトルで受け取る形式にしよう。

%評価領域、コイル、共に広大だが、どちらの計算をベクトル化した方が良いかゆっくり考えた方が良い。
%理想的には両方をベクトル化するのが望ましいが、コイル1e4エレメント、評価点6400点の場合、64000k点の情報を保持する必要がある。
%評価点の数を減らして、コイルをベクトル化しよう。そっちの方が、最適化という点では良いだろう。
%あと、評価面をベクトル化したバージョンも作りましょう。面で評価したい場合に役に立つ。というか、まず評価面をベクトル化した方法で計算しますか。


ap=zeros(FieldSize^2,1);
bp=zeros(FieldSize^2,1);
zpp=zeros(FieldSize^2,1);

i=1;
for a=1:FieldSize
    for b=1:FieldSize
        
        ap(i)=a*1e-3;
        bp(i)=b*1e-3;
        zpp(i)=zp*1e-3;
        i=i+1;
        
    end
end


for i=1:cnt(3)
    j=1;
    
    while(cnt(1)>j)
       
        if A(j+1,1,i)==0 && A(j+1,1,i)==0 
            break
        end
        
        %plot3(A(j:j+1,1,i),A(j:j+1,2,i),[zc,zc],"r");
        
        
        x=A(j,1,i)*1e-3;
        y=A(j,2,i)*1e-3;
        xd=(A(j,1,i)-A(j+1,1,i))*1e-3;
        yd=(A(j,2,i)-A(j+1,2,i))*1e-3;
        
        %ベクトルで演算できるように改良しましょう。

             
            %Pp=[ap*1e-3,bp*1e-3,zpp*1e-3];
            
            %Pc=[x,y,zc*1e-3];

            dS=[xd,yd,zc*1e-3];

            %r=Pp-Pc;
            r=[ap-x,bp-y,zpp-zc*1e-3];
            %ノルムの3乗ベクトルを作った後に、cross(dS,r)を作って、dHとして、Fieldにぼんぼん加算していく方針でやりましょう。
            r_norm3=sum(r.^2,2);
            r_norm3=4*pi*r_norm3.^(1.5);
            
            
            %dH=I*cross(dS,r)/(4*pi*norm(r)^3);
            
            dH=I*(xd.*r(:,2)-yd.*r(:,1))./r_norm3;
            %dH=0;
            
            Field_vec=Field_vec+dH;
     
        %plot3(A(j:j+1,1,i),A(j:j+1,2,i),[-50,-50],"r");
        

        %hold on;
        j=j+1;
    end
    i
end
% 

for i=1:FieldSize
    
    Field(:,i)=Field_vec(1+FieldSize*(i-1):FieldSize*i);
    
end

f=0;
end

