%Coil data gene

%��`�f�ʃR�C���f�[�^�쐬
function Coil=CoilGene(rp,rb,ac,ap,Num_lay,lpn)


%zc_w=0;             %�ŉ��w�R�C���̍Ŏ�O���W
%rp=20*1e-3;         %�v���W�F�N�^�C���O�a
%ac=1*1e-3;          %�ŉ��w�R�C������������
%rb=1*1e-3;          %�o������
rc=rp+rb+ac*0.5+ap*0.5;        %�ŉ��w�R�C�����a
% Total_turn=0;

%Num_lay=4;          %Number of layer
%lpn=10;             %Number of lay per num

%�R�C���ɕK�v�ȃf�[�^���o����
%�R�C����������Nc�Ƃ���ƁAi��Nc�܂łӂ������ɁA���ׂẴR�C���ɃA�N�Z�X�ł���K�v������B
%�܂�A�R�C���̎w��́A�P��̃A�h���X�ň����o���K�v������
%�K�v�ȃf�[�^�́A���a�f�[�^�A���C���a�Az���W�ł���B
%z���W�́Azc_w����̑��Έʒu��������Ă���΂悢�B
%���a�́A���C���a���l�����āA�����Ȃ��悤�B
%�R�C���f�[�^���v���b�g���āA�������Ȃ������ׂĂ݂悤�B

Coil=zeros(Num_lay*lpn,3);

for i=0:Num_lay-1

    for j=0:lpn-1
    
       Coil(i*lpn+j+1,1)=rc+ac*i;
       Coil(i*lpn+j+1,2)=ac*j;      %���̒l��zc_w�����Z����ƁA���[���h���W�ɂȂ�B
       Coil(i*lpn+j+1,3)=ac;

    end
end

end

