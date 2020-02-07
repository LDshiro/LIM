x_num=512;
y_num=512;
Amp=1e-6;
%a=107;      %zc [mm]
a=65;      %zc [mm]


F=transpose(C);
point=1;
line_size=size(F(1:end,1));

max_size=0;

cnt=1;

while((point-1)~=line_size(1))

    Next=F(point,2);
    point=point+Next+1;
    
    cnt=cnt+1;
    
    if max_size<Next
        
        max_size=Next;
        
    end
    
end

A=zeros(1500,2,100);
B=zeros(1500,2,100);
point=1;

for i=1:cnt-1

    Next=F(point,2);
    Value=F(point,1);
    
    hh=F(point+Next,:)-F(point+1,:);
    h=F(point+1,2);
    
    if vecnorm(hh)<2
        
        if Value > 0
            
            %if h>256/2
                A(1:Next,:,i)=F(point+1:(point+Next),:);
            %end
           
            %もう一つ格納するようのベクトルとして、Bベクトルを作って、始点の位置（hh）が,コイル中心の上にあるか下にあるかで、切り分けよう
        else
            %if h<256/2
                B(1:Next,:,i)=F(point+1:(point+Next),:);
            %end
            
        end
        
    end
    
    
    point=point+Next+1;

    
    if point > line_size(1)
        break
    end
    
end

   DD=BioSava_gradient2(A,a,0,256,1)+BioSava_gradient2(A,-a,0,256,1);
%   half=0
   CC=BioSava_gradient2(B,a,0,256,1)+BioSava_gradient2(B,-a,0,256,1);

     %DD=BioSava_gradientZ_mex(A,150,0,256,1)-BioSava_gradientZ_mex(A,-150,0,256,1);
     %CC=BioSava_gradientZ_mex(B,150,0,256,1)-BioSava_gradientZ_mex(B,-150,0,256,1);

%ここではAのコイルをプロットしよう

%max_size=50;
fa=0;


for i=1:cnt
    j=1;
    
    while(max_size>j)
       
        if A(j+1,1,i)==0 && A(j+1,1,i)==0 
            break
        end
        
        plot3(A(j:j+1,1,i),A(j:j+1,2,i),[50,50],"r");
        plot3(A(j:j+1,1,i),A(j:j+1,2,i),[-50,-50],"r");
        
        %line(A(j+1,1,i),A(j+1,2,i));
        hold on;
        j=j+1;
    end
    
end
% 

for i=1:cnt
    j=1;
    
    while(max_size>j)
       
        if B(j+1,1,i)==0 && B(j+1,1,i)==0 
            break
        end
        
        plot3(B(j:j+1,1,i),B(j:j+1,2,i),[50,50],"b");
        plot3(B(j:j+1,1,i),B(j:j+1,2,i),[-50,-50],"b");
        
        %line(A(j+1,1,i),A(j+1,2,i));
        hold on;
        j=j+1;
    end
    
end

%subplot(2,1,1);
K=DD+CC;
contourf(DD+CC,30);
 zlim([-70,70])

% subplot(2,1,2);
hold on;

for i=1:10
    plot(K(1:end,i*5+100));
end


% 
% [X,Y,Z] = meshgrid(-100:4:100);
% V = X.*exp(-X.^2-Y.^2-Z.^2);
% 
% xslice = [-100,0,100];   
% yslice = [];
% zslice = 0;
% slice(X,Y,Z,V,xslice,yslice,zslice)
