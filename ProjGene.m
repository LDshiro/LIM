%Projectle data gene

function Proj=ProjGene(rp,ap,Num_lay,lpn)
%Proj(�ŊO���a�A����\�A�w���A�w�����芪����)
%��`�f�ʃv���W�F�N�^�C���f�[�^�쐬

%rp=20*1e-3;         %�X���[�u�ŊO�a
%ap=1*1e-3;          %�X���[�u����\
%Num_lay=3;          %Number of layer
%lpn=10;             %Number of lay per num

%�R�C���ɕK�v�ȃf�[�^���o�����Ɠ���
%�R�C����������Nc�Ƃ���ƁAi��Nc�܂łӂ������ɁA���ׂẴR�C���ɃA�N�Z�X�ł���K�v������B
%�܂�A�R�C���̎w��́A�P��̃A�h���X�ň����o���K�v������
%�K�v�ȃf�[�^�́A���a�f�[�^�A���C���a�Az���W�ł���B
%z���W�́Azc_w����̑��Έʒu��������Ă���΂悢�B
%���a�́A���C���a���l�����āA�����Ȃ��悤�B
%�R�C���f�[�^���v���b�g���āA�������Ȃ������ׂĂ݂悤�B

Nc=Num_lay*lpn;

Proj=zeros(Nc,3);

for i=0:Num_lay-1

    for j=0:lpn-1
    
       Proj(i*lpn+j+1,1)=rp-ap*i;   %�ŊO�a����ɁA�����ɐϑw���銴��
       Proj(i*lpn+j+1,2)=ap*j;      %���̒l��zp_w�����Z����ƁA���[���h���W�ɂȂ�B
       Proj(i*lpn+j+1,3)=ap;        %���C�������i���a
%        Proj(i*lpn+j+1,4)=R_Al(Proj(i*lpn+j+1,1)^2*pi,ap,ap,20);         %R_Al(length,H,W,temp)
%        Proj(i*lpn+j+1,5)=SelfInd(Proj(i*lpn+j+1,1),Proj(i*lpn+j+1,3)*0.5);  %SelfInd(���a,���C�����a);
    end
end


end