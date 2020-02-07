%gradienr coil mapper

function Field = BioSava_gradient2(A,zc,zp,FieldSize,I)

Field=zeros(FieldSize,FieldSize);
Field_vec=zeros(FieldSize*FieldSize,1);

cnt=size(A);

%���ǂ̕����Ƃ��ẮAA��A���̃x�N�g���ɕϊ���������Ŏ��g�����B
%��?�ł��AA�͑O��̍��������x�N�g���ł�����̂ł�?
%�����ƁA

% x=A(j,1,i)*1e-3;
% y=A(j,2,i)*1e-3;
% xd=(A(j,1,i)-A(j+1,1,i))*1e-3;
% yd=(A(j,2,i)-A(j+1,2,i))*1e-3;

%����i�ɑ΂��āA���̏���������4��̃x�N�g����p�ӂ��Ă�����΁A�ʓ|�Ȍv�Z����C�ɏ􍞂߂�B

%�]���̈���x�N�g���Ŏ󂯎��`���ɂ��悤�B

%�]���̈�A�R�C���A���ɍL�傾���A�ǂ���̌v�Z���x�N�g�������������ǂ����������l���������ǂ��B
%���z�I�ɂ͗������x�N�g��������̂��]�܂������A�R�C��1e4�G�������g�A�]���_6400�_�̏ꍇ�A64000k�_�̏���ێ�����K�v������B
%�]���_�̐������炵�āA�R�C�����x�N�g�������悤�B�������̕����A�œK���Ƃ����_�ł͗ǂ����낤�B
%���ƁA�]���ʂ��x�N�g���������o�[�W���������܂��傤�B�ʂŕ]���������ꍇ�ɖ��ɗ��B�Ƃ������A�܂��]���ʂ��x�N�g�����������@�Ōv�Z���܂����B


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
        
        %�x�N�g���ŉ��Z�ł���悤�ɉ��ǂ��܂��傤�B

             
            %Pp=[ap*1e-3,bp*1e-3,zpp*1e-3];
            
            %Pc=[x,y,zc*1e-3];

            dS=[xd,yd,zc*1e-3];

            %r=Pp-Pc;
            r=[ap-x,bp-y,zpp-zc*1e-3];
            %�m������3��x�N�g�����������ɁAcross(dS,r)������āAdH�Ƃ��āAField�ɂڂ�ڂ���Z���Ă������j�ł��܂��傤�B
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

