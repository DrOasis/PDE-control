clear;
clc;
load RDNN_DSS_5R_2D_v1 y_y N_t

Z=[];
for iz=0:0.3:3
    for iq=0:0.3:3
         Z=[Z;iz];
    end
end
ZZ=[];

for it=0:20
     ZZ=[ZZ;Z];
end

Q=[];
for iz=0:0.3:3
    for iq=0:0.3:3
         Q=[Q;iq];
    end
end
QQ=[];
for it=0:20
     QQ=[QQ;Q];
end

TT=[];
for it=0:0.5:10
    for iz=0:0.3:3
        for iq=0:0.3:3
             TT=[TT;it];
        end
    end
end
m1=zeros(11,11,N_t);

for it=1:N_t
for iz=1:11
    for iq=1:11

      m1(iz,iq,it)=y_y((6*iz)-5,(6*iq)-5,it);
     end
end
end

YY=[];
for it=0:20
     if it==0                
        M1=m1(:,:,it+1);
     else
M1=m1(:,:,it*50);
    end
     
     YY=[YY;M1(:)];
end

save Data_writing_close ZZ QQ YY TT


