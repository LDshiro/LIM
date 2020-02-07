clear

zc_w=0;             %�ŉ��w�R�C���̍Ŏ�O���W
rp=20*1e-3;         %�v���W�F�N�^�C���O�a
ac=1*1e-3;          %�ŉ��w�R�C������������
ap=1*1e-3;          %�ŉ��w�R�C������������

rb=1*1e-3;          %�o������
rc=rp+rb+ac;        %�ŉ��w�R�C�����a
Num_layC=4;          %Number of layer
Num_lay=1;
lpn=35;             %Number of lay per num

lpnC=17;
z_start=0.0120;

%CoilGene�@�\�ɁA�R�C���p�����[�^�z����o�͂���@�\��t�����܂��傤�B

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


%%A�s��̒萔�̓���
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


%%����A�萔����ˁH���萔�ł��킱��
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


%�����I�ɏ����Ă�
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




