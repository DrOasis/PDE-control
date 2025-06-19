%%%RDNN_DSS_5R_2D
%%%2D space，5 robots

clear;
clc;

D=[0.01 0;
    0 0.01];
F=0.1*[1 1];

L=1;
lambda=1.5;
chi=0.1;

setlmis([])
P1=lmivar(1,[1 1]);
P2=lmivar(1,[1 1]);
Khat=lmivar(1,[1 1]);

fir=1;
lmiterm([fir 1 1 Khat],-1,L);
lmiterm([fir 1 1 Khat'],-L,1);
lmiterm([fir 1 1 P1],lambda,1);
lmiterm([fir 1 2 P1],1,1);
lmiterm([fir 1 2 P2],-1,1);
lmiterm([fir 1 2 Khat'],-L',1);
lmiterm([fir 1 3 P2],-1,F(1,1));
lmiterm([fir 1 3 Khat'],-L',F(1,1));
lmiterm([fir 1 4 P2],-1,F(1,2));
lmiterm([fir 1 4 Khat'],-L',F(1,2));
lmiterm([fir 1 6 Khat],-1,L);

lmiterm([fir 2 2 P2],-2,1);
lmiterm([fir 2 3 P2],-2,F(1,1));
lmiterm([fir 2 4 P2],-2,F(1,2));
lmiterm([fir 2 5 P2],1,D(1,1)+D(2,2));
lmiterm([fir 2 6 Khat],-1,L);

lmiterm([fir 3 3 P2],-2*F(1,1),F(1,1));
lmiterm([fir 3 4 P2],-F(1,1),F(1,2));
lmiterm([fir 4 4 P2],-2*F(1,2),F(1,2));

lmiterm([fir 3 3 P2],lambda,D(1,1));
lmiterm([fir 3 3 P2],-2,D(1,1));
lmiterm([fir 4 4 P2],-2,D(2,2));
lmiterm([fir 3 3 Khat],-2*D(1,1),L);
lmiterm([fir 4 4 Khat],-2*D(2,2),L);

lmiterm([fir 3 5 P2],2*F(1,1),D(1,1)+D(2,2));
lmiterm([fir 4 5 P2],2*F(1,2),D(1,1)+D(2,2));

lmiterm([fir 3 6 Khat],-F(1,1),L);
lmiterm([fir 4 6 Khat],-F(1,2),L);

lmiterm([fir 5 5 P2],-2,D(1,1)+D(2,2));
lmiterm([fir 5 6 Khat],D(1,1)+D(2,2),L);

lmiterm([fir 6 6 P2],-lambda*pi^2/(4*chi^2),D(2,2));

fir=fir+1;
lmiterm([fir 1 1 P1],-1,1);
fir=fir+1;
lmiterm([fir 1 1 P2],-1,1);

lmis=getlmis;
[tmin,xfeas]=feasp(lmis);
P1=dec2mat(lmis,xfeas,P1);
P2=dec2mat(lmis,xfeas,P2);
Khat=dec2mat(lmis,xfeas,Khat);
K_o=P2^-1*Khat;

a=0.2;
k_u=6;
k_v=4;

t0=0;
tf=10;
% tf=5;
t_sample=0.01;
N_t=round((tf-t0)/t_sample)+1;

z0=0;
z1=3;
z_sample=0.05;
N_z=round((z1-z0)/z_sample)+1;

q0=0;
q1=3;
q_sample=0.05;
N_q=round((q1-q0)/q_sample)+1;

ia=round(a/z_sample);

ttt=zeros(N_t,1);
zzz=zeros(N_z,1);
qqq=zeros(N_q,1);
for it=1:N_t
    ttt(it)=(it-1)*t_sample;
end
for iz=1:N_z
    zzz(iz)=(iz-1)*z_sample;
end
for iq=1:N_q
    qqq(iq)=(iq-1)*q_sample;
end

y=zeros(N_t,N_z,N_q);
dy_z=zeros(N_t,N_z,N_q);
ddy_z=zeros(N_t,N_z,N_q);
dy_q=zeros(N_t,N_z,N_q);
ddy_q=zeros(N_t,N_z,N_q);
uk=zeros(N_t,N_z,N_q);
v1=zeros(N_t,2);
v2=zeros(N_t,2);
v3=zeros(N_t,2);
v4=zeros(N_t,2);
v5=zeros(N_t,2);
s1=zeros(N_t,2);
s2=zeros(N_t,2);
s3=zeros(N_t,2);
s4=zeros(N_t,2);
s5=zeros(N_t,2);
y1=zeros(N_t,N_z,N_q);
error=zeros(N_t,N_z,N_q);
dy1_z=zeros(N_t,N_z,N_q);
ddy1_z=zeros(N_t,N_z,N_q);
dy1_q=zeros(N_t,N_z,N_q);
ddy1_q=zeros(N_t,N_z,N_q);
u_o=zeros(N_t,N_z,N_q);

N=5;
uk_hat=zeros(N_t+N,N_z,N_q);
s1_hat=zeros(N_t+N,2);
s2_hat=zeros(N_t+N,2);
s3_hat=zeros(N_t+N,2);
s4_hat=zeros(N_t+N,2);
s5_hat=zeros(N_t+N,2);
v1_hat=zeros(N_t+N,2);
v2_hat=zeros(N_t+N,2);
v3_hat=zeros(N_t+N,2);
v4_hat=zeros(N_t+N,2);
v5_hat=zeros(N_t+N,2);

s10=[1.8 0.5];
is10=[round(s10(1,1)/z_sample)+1 round(s10(1,2)/q_sample)+1];
s1(1,1)=s10(1,1);
s1(1,2)=s10(1,2);

s20=[2.5 1.8];
is20=[round(s20(1,1)/z_sample)+1 round(s20(1,2)/q_sample)+1];
s2(1,1)=s20(1,1);
s2(1,2)=s20(1,2);

s30=[1.2 2.5];
is30=[round(s30(1,1)/z_sample)+1 round(s30(1,2)/q_sample)+1];
s3(1,1)=s30(1,1);
s3(1,2)=s30(1,2);

s40=[0.5 1];
is40=[round(s40(1,1)/z_sample)+1 round(s40(1,2)/q_sample)+1];
s4(1,1)=s40(1,1);
s4(1,2)=s40(1,2);

s50=[1.5 1];
is50=[round(s50(1,1)/z_sample)+1 round(s50(1,2)/q_sample)+1];
s5(1,1)=s50(1,1);
s5(1,2)=s50(1,2);

for iz=1:N_z
    for iq=1:N_q

    y(1,iz,iq)=1*sin(pi*(iz-1)*z_sample)*sin(pi*(iq-1)*q_sample);
    dy_z(1,iz,iq)=1*pi*cos(pi*(iz-1)*z_sample)*sin(pi*(iq-1)*q_sample);
    ddy_z(1,iz,iq)=-1*pi*pi*sin(pi*(iz-1)*z_sample)*sin(pi*(iq-1)*q_sample);  
    dy_q(1,iz,iq)=1*pi*sin(pi*(iz-1)*z_sample)*cos(pi*(iq-1)*q_sample);
    ddy_q(1,iz,iq)=-1*pi*pi*sin(pi*(iz-1)*z_sample)*sin(pi*(iq-1)*q_sample);
        if y(1,iz,iq)<0
             y(1,iz,iq)=0;
        end

y1(1,iz,iq)=0;
    dy1_z(1,iz,iq)=0;
    ddy1_z(1,iz,iq)=0;
    dy1_q(1,iz,iq)=0;
    ddy1_q(1,iz,iq)=0;
    
error(1,iz,iq)= y(1,iz,iq)-y1(1,iz,iq);
if 1<=iq && iq<3; u_o(1,iz,iq)=K_o*L*error(1,iz,2);
elseif 3<=iq && iq<5; u_o(1,iz,iq)=K_o*L*error(1,iz,4);
elseif 5<=iq && iq<7; u_o(1,iz,iq)=K_o*L*error(1,iz,6);
elseif 7<=iq && iq<9; u_o(1,iz,iq)=K_o*L*error(1,iz,8);
elseif 9<=iq && iq<11; u_o(1,iz,iq)=K_o*L*error(1,iz,10);
elseif 11<=iq && iq<13; u_o(1,iz,iq)=K_o*L*error(1,iz,12);
elseif 13<=iq && iq<15; u_o(1,iz,iq)=K_o*L*error(1,iz,14);
elseif 15<=iq && iq<17; u_o(1,iz,iq)=K_o*L*error(1,iz,16);
elseif 17<=iq && iq<19; u_o(1,iz,iq)=K_o*L*error(1,iz,18);
elseif 19<=iq && iq<21; u_o(1,iz,iq)=K_o*L*error(1,iz,20);
elseif 21<=iq && iq<23; u_o(1,iz,iq)=K_o*L*error(1,iz,22);
elseif 23<=iq && iq<25; u_o(1,iz,iq)=K_o*L*error(1,iz,24);
elseif 25<=iq && iq<27; u_o(1,iz,iq)=K_o*L*error(1,iz,26);
elseif 27<=iq && iq<29; u_o(1,iz,iq)=K_o*L*error(1,iz,28);
elseif 29<=iq && iq<31; u_o(1,iz,iq)=K_o*L*error(1,iz,30); 
elseif 31<=iq && iq<33; u_o(1,iz,iq)=K_o*L*error(1,iz,32);
elseif 33<=iq && iq<35; u_o(1,iz,iq)=K_o*L*error(1,iz,34);
elseif 35<=iq && iq<37; u_o(1,iz,iq)=K_o*L*error(1,iz,36);
elseif 37<=iq && iq<39; u_o(1,iz,iq)=K_o*L*error(1,iz,38);
elseif 39<=iq && iq<41; u_o(1,iz,iq)=K_o*L*error(1,iz,40);
elseif 41<=iq && iq<43; u_o(1,iz,iq)=K_o*L*error(1,iz,42);
elseif 43<=iq && iq<45; u_o(1,iz,iq)=K_o*L*error(1,iz,44);
elseif 45<=iq && iq<47; u_o(1,iz,iq)=K_o*L*error(1,iz,46);
elseif 47<=iq && iq<49; u_o(1,iz,iq)=K_o*L*error(1,iz,48);
elseif 49<=iq && iq<51; u_o(1,iz,iq)=K_o*L*error(1,iz,50);
elseif 51<=iq && iq<53; u_o(1,iz,iq)=K_o*L*error(1,iz,52);
elseif 53<=iq && iq<55; u_o(1,iz,iq)=K_o*L*error(1,iz,54);
elseif 55<=iq && iq<57; u_o(1,iz,iq)=K_o*L*error(1,iz,56);
elseif 57<=iq && iq<59; u_o(1,iz,iq)=K_o*L*error(1,iz,58);
elseif 59<=iq && iq<=61; u_o(1,iz,iq)=K_o*L*error(1,iz,60);
end
    end
end

for it=1:N_t
    for iq=1:N_q
     y(it,1,iq)=0;
     y(it,N_z,iq)=0;
     y1(it,1,iq)=0;
     y1(it,N_z,iq)=0;
    end
    for iz=1:N_z
     y(it,iz,1)=0;
     y(it,iz,N_q)=0;
     y1(it,iz,1)=0;
     y1(it,iz,N_q)=0;
    end
end



A10=2*a*a*(squeeze(y(1,is10(1,1)-ia,is10(1,2)-ia))+squeeze(y(1,is10(1,1)+ia,is10(1,2)-ia)))/2 ...
                +(squeeze(y(1,is10(1,1)-ia,is10(1,2)+ia))+squeeze(y(1,is10(1,1)+ia,is10(1,2)+ia)))/2;
A20=2*a*a*(squeeze(y(1,is20(1,1)-ia,is20(1,2)-ia))+squeeze(y(1,is20(1,1)+ia,is20(1,2)-ia)))/2 ...
                +(squeeze(y(1,is20(1,1)-ia,is20(1,2)+ia))+squeeze(y(1,is20(1,1)+ia,is20(1,2)+ia)))/2;
A30=2*a*a*(squeeze(y(1,is30(1,1)-ia,is30(1,2)-ia))+squeeze(y(1,is30(1,1)+ia,is30(1,2)-ia)))/2 ...
                +(squeeze(y(1,is30(1,1)-ia,is30(1,2)+ia))+squeeze(y(1,is30(1,1)+ia,is30(1,2)+ia)))/2;
A40=2*a*a*(squeeze(y(1,is40(1,1)-ia,is40(1,2)-ia))+squeeze(y(1,is40(1,1)+ia,is40(1,2)-ia)))/2 ...
                +(squeeze(y(1,is40(1,1)-ia,is40(1,2)+ia))+squeeze(y(1,is40(1,1)+ia,is40(1,2)+ia)))/2;
A50=2*a*a*(squeeze(y(1,is50(1,1)-ia,is50(1,2)-ia))+squeeze(y(1,is50(1,1)+ia,is50(1,2)-ia)))/2 ...
                +(squeeze(y(1,is50(1,1)-ia,is50(1,2)+ia))+squeeze(y(1,is50(1,1)+ia,is50(1,2)+ia)))/2;

ur1(1)=-k_u*A10;
ur2(1)=-k_u*A20;
ur3(1)=-k_u*A30;
ur4(1)=-k_u*A40;
ur5(1)=-k_u*A50;

v1(1,1)=-ur1(1)*k_v*(squeeze(y(1,is10(1,1)+ia,is10(1,2)))-squeeze(y(1,is10(1,1)-ia,is10(1,2))))/(2*a);
v1(1,2)=-ur1(1)*k_v*(squeeze(y(1,is10(1,1),is10(1,2)+ia))-squeeze(y(1,is10(1,1),is10(1,2)-ia)))/(2*a);
v2(1,1)=-ur2(1)*k_v*(squeeze(y(1,is20(1,1)+ia,is20(1,2)))-squeeze(y(1,is20(1,1)-ia,is20(1,2))))/(2*a);
v2(1,2)=-ur2(1)*k_v*(squeeze(y(1,is20(1,1),is20(1,2)+ia))-squeeze(y(1,is20(1,1),is20(1,2)-ia)))/(2*a);
v3(1,1)=-ur3(1)*k_v*(squeeze(y(1,is30(1,1)+ia,is30(1,2)))-squeeze(y(1,is30(1,1)-ia,is30(1,2))))/(2*a);
v3(1,2)=-ur3(1)*k_v*(squeeze(y(1,is30(1,1),is30(1,2)+ia))-squeeze(y(1,is30(1,1),is30(1,2)-ia)))/(2*a);
v4(1,1)=-ur4(1)*k_v*(squeeze(y(1,is40(1,1)+ia,is40(1,2)))-squeeze(y(1,is40(1,1)-ia,is40(1,2))))/(2*a);
v4(1,2)=-ur4(1)*k_v*(squeeze(y(1,is40(1,1),is40(1,2)+ia))-squeeze(y(1,is40(1,1),is40(1,2)-ia)))/(2*a);
v5(1,1)=-ur5(1)*k_v*(squeeze(y(1,is50(1,1)+ia,is50(1,2)))-squeeze(y(1,is50(1,1)-ia,is50(1,2))))/(2*a);
v5(1,2)=-ur5(1)*k_v*(squeeze(y(1,is50(1,1),is50(1,2)+ia))-squeeze(y(1,is50(1,1),is50(1,2)-ia)))/(2*a);

q=1;
Q_cell=cell(N,N);
for i=1:1:N
    for j=1:1:N
        if i==j
            Q_cell{i,j}=q;
        else
            Q_cell{i,j}=zeros(1);
        end
    end
end
Q=1*cell2mat(Q_cell);
gamma=1;
D_n=1.5;
D_hat=kron(Q,D_n);
Qpinv=pinv(Q);

xN=zeros(N_t,N_z,N_q,N);
dxN_z=zeros(N_t,N_z,N_q,N);
dxN_q=zeros(N_t,N_z,N_q,N);
ddxN_z=zeros(N_t,N_z,N_q,N);
ddxN_q=zeros(N_t,N_z,N_q,N);
ddxN=zeros(N_t,N_z,N_q,N);
dddxN_z=zeros(N_t,N_z,N_q,N);
dddxN_q=zeros(N_t,N_z,N_q,N);
ddddxN_z=zeros(N_t,N_z,N_q,N);
ddddxN_q=zeros(N_t,N_z,N_q,N);
dddxN_t=zeros(N_t,N_z,N_q,N);
dU_t=zeros(N_t,N_z,N_q,N);
DDxN=zeros(N_t,N_z,N_q,N);
DxN_z=zeros(N_t,N_z,N_q,N);
DxN_q=zeros(N_t,N_z,N_q,N);
dxN_zq=zeros(N_t,N_z,N_q,N);
ddxN_zq_t=zeros(N_t,N_z,N_q,N);
dcF_z=zeros(N_t,N_z,N_q,N);
dcF_q=zeros(N_t,N_z,N_q,N);
ddcF_z=zeros(N_t,N_z,N_q,N);
ddcF_q=zeros(N_t,N_z,N_q,N);

U=zeros(N_t,N_z,N_q,N);
y1_hat=zeros(N_t+N,N_z,N_q);
dy1_hat_z=zeros(N_t+N,N_z,N_q);
ddy1_hat_z=zeros(N_t+N,N_z,N_q);
dy1_hat_q=zeros(N_t+N,N_z,N_q);
ddy1_hat_q=zeros(N_t+N,N_z,N_q);

A1_hat=zeros(N_t+N,1);
A2_hat=zeros(N_t+N,1);
A3_hat=zeros(N_t+N,1);
A4_hat=zeros(N_t+N,1);
A5_hat=zeros(N_t+N,1);
ur1_hat=zeros(N_t+N,1);
ur2_hat=zeros(N_t+N,1);
ur3_hat=zeros(N_t+N,1);
ur4_hat=zeros(N_t+N,1);
ur5_hat=zeros(N_t+N,1);
A1=zeros(N_t,1);
A2=zeros(N_t,1);
A3=zeros(N_t,1);
A4=zeros(N_t,1);
A5=zeros(N_t,1);
ur1=zeros(N_t,1);
ur2=zeros(N_t,1);
ur3=zeros(N_t,1);
ur4=zeros(N_t,1);
ur5=zeros(N_t,1);
for it=1:N_t-1
    it
    
    is1=[round(s1(it,1)/z_sample) round(s1(it,2)/q_sample)];
    is2=[round(s2(it,1)/z_sample) round(s2(it,2)/q_sample)];
    is3=[round(s3(it,1)/z_sample) round(s3(it,2)/q_sample)];
    is4=[round(s4(it,1)/z_sample) round(s4(it,2)/q_sample)];
    is5=[round(s5(it,1)/z_sample) round(s5(it,2)/q_sample)];
    
    for iz=2:N_z-1  
        for iq=2:N_q-1
         y(it+1,iz,iq)=squeeze(y(it,iz,iq))+t_sample*(squeeze(D(1,1)*ddy_z(it,iz,iq)+D(2,2)*ddy_q(it,iz,iq))-squeeze(F(1,1)*dy_z(it,iz,iq)+F(1,2)*dy_q(it,iz,iq))+squeeze(uk(it,iz,iq)));

%          y(it+1,iz,iq)=squeeze(y(it,iz,iq))+t_sample*(squeeze(D(1,1)*ddy_z(it,iz,iq)+D(2,2)*ddy_q(it,iz,iq))-squeeze(F(1,1)*dy_z(it,iz,iq)+F(1,2)*dy_q(it,iz,iq)));   

         y1(it+1,iz,iq)=squeeze(y1(it,iz,iq))+t_sample*(squeeze(D(1,1)*ddy1_z(it,iz,iq)+D(2,2)*ddy1_q(it,iz,iq))-squeeze(F(1,1)*dy1_z(it,iz,iq)...
                      +F(1,2)*dy1_q(it,iz,iq))+squeeze(u_o(it,iz,iq)));
                  
         if y(it+1,iz,iq)<0
             y(it+1,iz,iq)=0;
         end

         error(it+1,iz,iq)=y(it+1,iz,iq)-y1(it+1,iz,iq);
        
        end
    end
    
    
    for iz=2:N_z
        for iq=2:N_q
         dy_z(it+1,iz,iq)=(squeeze(y(it+1,iz,iq))-squeeze(y(it+1,iz-1,iq)))./z_sample;
         dy_q(it+1,iz,iq)=(squeeze(y(it+1,iz,iq))-squeeze(y(it+1,iz,iq-1)))./q_sample;
         
         dy1_z(it+1,iz,iq)=(squeeze(y1(it+1,iz,iq))-squeeze(y1(it+1,iz-1,iq)))./z_sample;
         dy1_q(it+1,iz,iq)=(squeeze(y1(it+1,iz,iq))-squeeze(y1(it+1,iz,iq-1)))./q_sample;
        end
    end
     for iz=1:N_z-1
         for iq=1:N_q-1
         ddy_z(it+1,iz,iq)=(squeeze(dy_z(it+1,iz+1,iq))-squeeze(dy_z(it+1,iz,iq)))./z_sample;
         ddy_q(it+1,iz,iq)=(squeeze(dy_q(it+1,iz,iq+1))-squeeze(dy_q(it+1,iz,iq)))./q_sample;
         
         ddy1_z(it+1,iz,iq)=(squeeze(dy1_z(it+1,iz+1,iq))-squeeze(dy1_z(it+1,iz,iq)))./z_sample;
         ddy1_q(it+1,iz,iq)=(squeeze(dy1_q(it+1,iz,iq+1))-squeeze(dy1_q(it+1,iz,iq)))./q_sample;
         
if 1<=iq && iq<3; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,2);
elseif 3<=iq && iq<5; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,4);
elseif 5<=iq && iq<7; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,6);
elseif 7<=iq && iq<9; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,8);
elseif 9<=iq && iq<11; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,10);
elseif 11<=iq && iq<13; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,12);
elseif 13<=iq && iq<15; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,14);
elseif 15<=iq && iq<17; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,16);
elseif 17<=iq && iq<19; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,18);
elseif 19<=iq && iq<21; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,20);
elseif 21<=iq && iq<23; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,22);
elseif 23<=iq && iq<25; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,24);
elseif 25<=iq && iq<27; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,26);
elseif 27<=iq && iq<29; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,28);
elseif 29<=iq && iq<31; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,30); 
elseif 31<=iq && iq<33; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,32);
elseif 33<=iq && iq<35; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,34);
elseif 35<=iq && iq<37; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,36);
elseif 37<=iq && iq<39; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,38);
elseif 39<=iq && iq<41; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,40);
elseif 41<=iq && iq<43; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,42);
elseif 43<=iq && iq<45; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,44);
elseif 45<=iq && iq<47; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,46);
elseif 47<=iq && iq<49; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,48);
elseif 49<=iq && iq<51; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,50);
elseif 51<=iq && iq<53; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,52);
elseif 53<=iq && iq<55; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,54);
elseif 55<=iq && iq<57; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,56);
elseif 57<=iq && iq<59; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,58);
elseif 59<=iq && iq<=61; u_o(it+1,iz,iq)=K_o*L*error(it+1,iz,60);
end

y1_hat(it,iz,iq)=y1(it,iz,iq);
uk_hat(it,iz,iq)=uk(it,iz,iq);
ddy1_hat_z(it,iz,iq)=ddy1_z(it,iz,iq);
ddy1_hat_q(it,iz,iq)=ddy1_q(it,iz,iq);
dy1_hat_z(it,iz,iq)=dy1_z(it,iz,iq); 
dy1_hat_q(it,iz,iq)=dy1_q(it,iz,iq);
         end
     end

          s1_hat(it,1)=s1(it,1);
          s1_hat(it,2)=s1(it,2);
          s2_hat(it,1)=s2(it,1);
          s2_hat(it,2)=s2(it,2);
          s3_hat(it,1)=s3(it,1);
          s3_hat(it,2)=s3(it,2);
          s4_hat(it,1)=s4(it,1);
          s4_hat(it,2)=s4(it,2);
          s5_hat(it,1)=s5(it,1);
          s5_hat(it,2)=s5(it,2);
          
          v1_hat(it,1)=v1(it,1);
          v1_hat(it,2)=v1(it,2);
          v2_hat(it,1)=v2(it,1);
          v2_hat(it,2)=v2(it,2);
          v3_hat(it,1)=v3(it,1);
          v3_hat(it,2)=v3(it,2);
          v4_hat(it,1)=v4(it,1);
          v4_hat(it,2)=v4(it,2);   
          v5_hat(it,1)=v5(it,1);
          v5_hat(it,2)=v5(it,2);   

 for itt=0:N
  
    s1_hat(it+itt+1,1)=s1_hat(it+itt,1)+t_sample*v1_hat(it+itt,1);
    s1_hat(it+itt+1,2)=s1_hat(it+itt,2)+t_sample*v1_hat(it+itt,2); 
    s2_hat(it+itt+1,1)=s2_hat(it+itt,1)+t_sample*v2_hat(it+itt,1);
    s2_hat(it+itt+1,2)=s2_hat(it+itt,2)+t_sample*v2_hat(it+itt,2);
    s3_hat(it+itt+1,1)=s3_hat(it+itt,1)+t_sample*v3_hat(it+itt,1);
    s3_hat(it+itt+1,2)=s3_hat(it+itt,2)+t_sample*v3_hat(it+itt,2);
    s4_hat(it+itt+1,1)=s4_hat(it+itt,1)+t_sample*v4_hat(it+itt,1);
    s4_hat(it+itt+1,2)=s4_hat(it+itt,2)+t_sample*v4_hat(it+itt,2);
    s5_hat(it+itt+1,1)=s5_hat(it+itt,1)+t_sample*v5_hat(it+itt,1);
    s5_hat(it+itt+1,2)=s5_hat(it+itt,2)+t_sample*v5_hat(it+itt,2);

    is1_hat=[round(s1_hat(it+itt+1,1)/z_sample) round(s1_hat(it+itt+1,2)/q_sample)];
    is2_hat=[round(s2_hat(it+itt+1,1)/z_sample) round(s2_hat(it+itt+1,2)/q_sample)];
    is3_hat=[round(s3_hat(it+itt+1,1)/z_sample) round(s3_hat(it+itt+1,2)/q_sample)];
    is4_hat=[round(s4_hat(it+itt+1,1)/z_sample) round(s4_hat(it+itt+1,2)/q_sample)];
    is5_hat=[round(s5_hat(it+itt+1,1)/z_sample) round(s5_hat(it+itt+1,2)/q_sample)];


if is1_hat(1,1)<=ia; is1_hat(1,1)=ia+1; 
end
if is1_hat(1,2)<=ia; is1_hat(1,2)=ia+1; 
end
if is2_hat(1,1)<=ia; is2_hat(1,1)=ia+1; 
end
if is2_hat(1,2)<=ia; is2_hat(1,2)=ia+1; 
end
if is3_hat(1,1)<=ia; is3_hat(1,1)=ia+1; 
end
if is3_hat(1,2)<=ia; is3_hat(1,2)=ia+1; 
end
if is4_hat(1,1)<=ia; is4_hat(1,1)=ia+1; 
end
if is4_hat(1,2)<=ia; is4_hat(1,2)=ia+1; 
end
if is5_hat(1,1)<=ia; is5_hat(1,1)=ia+1; 
end
if is5_hat(1,2)<=ia; is5_hat(1,2)=ia+1; 
end
if is1_hat(1,1)>=N_z-ia; is1_hat(1,1)=N_z-ia;
end
if is1_hat(1,2)>=N_q-ia; is1_hat(1,2)=N_q-ia;
end
if is2_hat(1,1)>=N_z-ia; is2_hat(1,1)=N_z-ia;
end
if is2_hat(1,2)>=N_q-ia; is2_hat(1,2)=N_q-ia;
end
if is3_hat(1,1)>=N_z-ia; is3_hat(1,1)=N_z-ia;
end
if is3_hat(1,2)>=N_q-ia; is3_hat(1,2)=N_q-ia;
end
if is4_hat(1,1)>=N_z-ia; is4_hat(1,1)=N_z-ia;
end
if is4_hat(1,2)>=N_q-ia; is4_hat(1,2)=N_q-ia;
end
if is5_hat(1,1)>=N_z-ia; is5_hat(1,1)=N_z-ia;
end
if is5_hat(1,2)>=N_q-ia; is5_hat(1,2)=N_q-ia;
end
    
     for iz=2:N_z-1  
        for iq=2:N_q-1  
         y1_hat(it+itt+1,iz,iq)=squeeze(y1_hat(it+itt,iz,iq))+t_sample*(squeeze(D(1,1)*ddy1_hat_z(it+itt,iz,iq)+D(2,2)*ddy1_hat_q(it+itt,iz,iq))...
             -squeeze(F(1,1)*dy1_hat_z(it+itt,iz,iq)+F(1,2)*dy1_hat_q(it+itt,iz,iq))...
             +squeeze(uk_hat(it+itt,iz,iq)));
        end
     end 


      A1_hat(it+itt+1)=2*a*a*(squeeze(y1_hat(it+itt+1,is1_hat(1,1)-ia,is1_hat(1,2)-ia))+squeeze(y1_hat(it+itt+1,is1_hat(1,1)+ia,is1_hat(1,2)-ia)))/2 ...
                +(squeeze(y1_hat(it+itt+1,is1_hat(1,1)-ia,is1_hat(1,2)+ia))+squeeze(y1_hat(it+itt+1,is1_hat(1,1)+ia,is1_hat(1,2)+ia)))/2;
      A2_hat(it+itt+1)=2*a*a*(squeeze(y1_hat(it+itt+1,is2_hat(1,1)-ia,is2_hat(1,2)-ia))+squeeze(y1_hat(it+itt+1,is2_hat(1,1)+ia,is2_hat(1,2)-ia)))/2 ...
                +(squeeze(y1_hat(it+itt+1,is2_hat(1,1)-ia,is2_hat(1,2)+ia))+squeeze(y1_hat(it+itt+1,is2_hat(1,1)+ia,is2_hat(1,2)+ia)))/2;
      A3_hat(it+itt+1)=2*a*a*(squeeze(y1_hat(it+itt+1,is3_hat(1,1)-ia,is3_hat(1,2)-ia))+squeeze(y1_hat(it+itt+1,is3_hat(1,1)+ia,is3_hat(1,2)-ia)))/2 ...
                +(squeeze(y1_hat(it+itt+1,is3_hat(1,1)-ia,is3_hat(1,2)+ia))+squeeze(y1_hat(it+itt+1,is3_hat(1,1)+ia,is3_hat(1,2)+ia)))/2;
      A4_hat(it+itt+1)=2*a*a*(squeeze(y1_hat(it+itt+1,is4_hat(1,1)-ia,is4_hat(1,2)-ia))+squeeze(y1_hat(it+itt+1,is4_hat(1,1)+ia,is4_hat(1,2)-ia)))/2 ...
                +(squeeze(y1_hat(it+itt+1,is4_hat(1,1)-ia,is4_hat(1,2)+ia))+squeeze(y1_hat(it+itt+1,is4_hat(1,1)+ia,is4_hat(1,2)+ia)))/2;
A5_hat(it+itt+1)=2*a*a*(squeeze(y1_hat(it+itt+1,is5_hat(1,1)-ia,is5_hat(1,2)-ia))+squeeze(y1_hat(it+itt+1,is5_hat(1,1)+ia,is5_hat(1,2)-ia)))/2 ...
                +(squeeze(y1_hat(it+itt+1,is5_hat(1,1)-ia,is5_hat(1,2)+ia))+squeeze(y1_hat(it+itt+1,is5_hat(1,1)+ia,is5_hat(1,2)+ia)))/2;

ur1_hat(it+itt+1)=-k_u*A1_hat(it+itt+1);
ur2_hat(it+itt+1)=-k_u*A2_hat(it+itt+1); 
ur3_hat(it+itt+1)=-k_u*A3_hat(it+itt+1);
ur4_hat(it+itt+1)=-k_u*A4_hat(it+itt+1);
ur5_hat(it+itt+1)=-k_u*A5_hat(it+itt+1);

v1_hat(it+itt+1,1)=-ur1_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is1_hat(1,1)+ia,is1_hat(1,2)))-squeeze(y1_hat(it+itt+1,is1_hat(1,1)-ia,is1_hat(1,2))))/(2*a);%%机器人1的z方向速度
v1_hat(it+itt+1,2)=-ur1_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is1_hat(1,1),is1_hat(1,2)+ia))-squeeze(y1_hat(it+itt+1,is1_hat(1,1),is1_hat(1,2)-ia)))/(2*a);%%机器人1的q方向速度
v2_hat(it+itt+1,1)=-ur2_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is2_hat(1,1)+ia,is2_hat(1,2)))-squeeze(y1_hat(it+itt+1,is2_hat(1,1)-ia,is2_hat(1,2))))/(2*a);%%机器人2的z方向速度
v2_hat(it+itt+1,2)=-ur2_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is2_hat(1,1),is2_hat(1,2)+ia))-squeeze(y1_hat(it+itt+1,is2_hat(1,1),is2_hat(1,2)-ia)))/(2*a);%%机器人2的q方向速度
v3_hat(it+itt+1,1)=-ur3_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is3_hat(1,1)+ia,is3_hat(1,2)))-squeeze(y1_hat(it+itt+1,is3_hat(1,1)-ia,is3_hat(1,2))))/(2*a);%%机器人3的z方向速度
v3_hat(it+itt+1,2)=-ur3_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is3_hat(1,1),is3_hat(1,2)+ia))-squeeze(y1_hat(it+itt+1,is3_hat(1,1),is3_hat(1,2)-ia)))/(2*a);%%机器人3的q方向速度
v4_hat(it+itt+1,1)=-ur4_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is4_hat(1,1)+ia,is4_hat(1,2)))-squeeze(y1_hat(it+itt+1,is4_hat(1,1)-ia,is4_hat(1,2))))/(2*a);%%机器人4的z方向速度
v4_hat(it+itt+1,2)=-ur4_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is4_hat(1,1),is4_hat(1,2)+ia))-squeeze(y1_hat(it+itt+1,is4_hat(1,1),is4_hat(1,2)-ia)))/(2*a);%%机器人4的q方向速度
v5_hat(it+itt+1,1)=-ur5_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is5_hat(1,1)+ia,is5_hat(1,2)))-squeeze(y1_hat(it+itt+1,is5_hat(1,1)-ia,is5_hat(1,2))))/(2*a);%%机器人4的z方向速度
v5_hat(it+itt+1,2)=-ur5_hat(it+itt+1)*k_v*(squeeze(y1_hat(it+itt+1,is5_hat(1,1),is5_hat(1,2)+ia))-squeeze(y1_hat(it+itt+1,is5_hat(1,1),is5_hat(1,2)-ia)))/(2*a);%%机器人4的q方向速度
 

    for iz=2:N_z
        for iq=2:N_q
         dy1_hat_z(it+itt+1,iz,iq)=(squeeze(y1_hat(it+itt+1,iz,iq))-squeeze(y1_hat(it+itt+1,iz-1,iq)))./z_sample;
         dy1_hat_q(it+itt+1,iz,iq)=(squeeze(y1_hat(it+itt+1,iz,iq))-squeeze(y1_hat(it+itt+1,iz,iq-1)))./q_sample;
        end
    end
     for iz=1:N_z-1
        for iq=1:N_q-1
         ddy1_hat_z(it+itt+1,iz,iq)=(squeeze(dy1_hat_z(it+itt+1,iz+1,iq))-squeeze(dy1_hat_z(it+itt+1,iz,iq)))./z_sample;
         ddy1_hat_q(it+itt+1,iz,iq)=(squeeze(dy1_hat_q(it+itt+1,iz,iq+1))-squeeze(dy1_hat_q(it+itt+1,iz,iq)))./q_sample;
         
      if is1_hat(1,1)-ia<=iz && iz<=is1_hat(1,1)+ia
        if is1_hat(1,2)-ia<=iq && iq<=is1_hat(1,2)+ia
        uk_hat(it+itt+1,iz,iq)=-k_u*A1_hat(it+itt+1);
        end
      end
       if is2_hat(1,1)-ia<=iz && iz<=is2_hat(1,1)+ia
        if is2_hat(1,2)-ia<=iq && iq<=is2_hat(1,2)+ia
        uk_hat(it+itt+1,iz,iq)=-k_u*A2_hat(it+itt+1);
        end
       end
       if is3_hat(1,1)-ia<=iz && iz<=is3_hat(1,1)+ia
        if is3_hat(1,2)-ia<=iq && iq<=is3_hat(1,2)+ia
        uk_hat(it+itt+1,iz,iq)=-k_u*A3_hat(it+itt+1);
        end
       end
       if is4_hat(1,1)-ia<=iz && iz<=is4_hat(1,1)+ia
        if is4_hat(1,2)-ia<=iq && iq<=is4_hat(1,2)+ia
        uk_hat(it+itt+1,iz,iq)=-k_u*A4_hat(it+itt+1);
        end
       end
      if is5_hat(1,1)-ia<=iz && iz<=is5_hat(1,1)+ia
        if is5_hat(1,2)-ia<=iq && iq<=is5_hat(1,2)+ia
        uk_hat(it+itt+1,iz,iq)=-k_u*A5_hat(it+itt+1);
        end
       end
     
       end
     end            

 end

      for iz=1:N_z-1
        for iq=1:N_q-1
        xN(it,iz,iq,1)=y1_hat(it+1,iz,iq);
        xN(it,iz,iq,2)=y1_hat(it+2,iz,iq);
        xN(it,iz,iq,3)=y1_hat(it+3,iz,iq);
        xN(it,iz,iq,4)=y1_hat(it+4,iz,iq);
        xN(it,iz,iq,5)=y1_hat(it+5,iz,iq);
        
        DDxN(it,iz,iq,1)=ddy1_hat_z(it+1,iz,iq)+ddy1_hat_q(it+1,iz,iq);
        DDxN(it,iz,iq,2)=ddy1_hat_z(it+2,iz,iq)+ddy1_hat_q(it+2,iz,iq);
        DDxN(it,iz,iq,3)=ddy1_hat_z(it+3,iz,iq)+ddy1_hat_q(it+3,iz,iq);
        DDxN(it,iz,iq,4)=ddy1_hat_z(it+4,iz,iq)+ddy1_hat_q(it+4,iz,iq);
        DDxN(it,iz,iq,5)=ddy1_hat_z(it+5,iz,iq)+ddy1_hat_q(it+5,iz,iq);
        
        DxN_z(it,iz,iq,1)=dy1_hat_z(it+1,iz,iq);
        DxN_z(it,iz,iq,2)=dy1_hat_z(it+2,iz,iq);
        DxN_z(it,iz,iq,3)=dy1_hat_z(it+3,iz,iq);
        DxN_z(it,iz,iq,4)=dy1_hat_z(it+4,iz,iq);
        DxN_z(it,iz,iq,5)=dy1_hat_z(it+5,iz,iq);
        DxN_q(it,iz,iq,1)=dy1_hat_q(it+1,iz,iq);
        DxN_q(it,iz,iq,2)=dy1_hat_q(it+2,iz,iq);
        DxN_q(it,iz,iq,3)=dy1_hat_q(it+3,iz,iq);
        DxN_q(it,iz,iq,4)=dy1_hat_q(it+4,iz,iq);
        DxN_q(it,iz,iq,5)=dy1_hat_q(it+5,iz,iq);

        ddxN_z(it,iz,iq,1)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+1,iz,iq));
        ddxN_z(it,iz,iq,2)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+2,iz,iq));
        ddxN_z(it,iz,iq,3)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+3,iz,iq));
        ddxN_z(it,iz,iq,4)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+4,iz,iq));
        ddxN_z(it,iz,iq,5)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+5,iz,iq));
        ddxN_z(it+1,iz,iq,1)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+2,iz,iq));
        ddxN_z(it+1,iz,iq,2)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+3,iz,iq));
        ddxN_z(it+1,iz,iq,3)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+4,iz,iq));
        ddxN_z(it+1,iz,iq,4)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+5,iz,iq));
        ddxN_z(it+1,iz,iq,5)=t_sample*D(1,1)*squeeze(ddy1_hat_z(it+6,iz,iq));
        
        ddxN_q(it,iz,iq,1)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+1,iz,iq));
        ddxN_q(it,iz,iq,2)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+2,iz,iq));
        ddxN_q(it,iz,iq,3)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+3,iz,iq));
        ddxN_q(it,iz,iq,4)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+4,iz,iq));
        ddxN_q(it,iz,iq,5)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+5,iz,iq));
        ddxN_q(it+1,iz,iq,1)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+2,iz,iq));
        ddxN_q(it+1,iz,iq,2)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+3,iz,iq));
        ddxN_q(it+1,iz,iq,3)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+4,iz,iq));
        ddxN_q(it+1,iz,iq,4)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+5,iz,iq));
        ddxN_q(it+1,iz,iq,5)=t_sample*D(2,2)*squeeze(ddy1_hat_q(it+6,iz,iq));
        
        ddxN(it,iz,iq,1)=ddxN_q(it,iz,iq,1)+ddxN_z(it,iz,iq,1);
        ddxN(it,iz,iq,2)=ddxN_q(it,iz,iq,2)+ddxN_z(it,iz,iq,2);
        ddxN(it,iz,iq,3)=ddxN_q(it,iz,iq,3)+ddxN_z(it,iz,iq,3);
        ddxN(it,iz,iq,4)=ddxN_q(it,iz,iq,4)+ddxN_z(it,iz,iq,4);
        ddxN(it,iz,iq,5)=ddxN_q(it,iz,iq,5)+ddxN_z(it,iz,iq,5);
        ddxN(it+1,iz,iq,1)=ddxN_q(it+1,iz,iq,1)+ddxN_z(it+1,iz,iq,1);
        ddxN(it+1,iz,iq,2)=ddxN_q(it+1,iz,iq,2)+ddxN_z(it+1,iz,iq,2);
        ddxN(it+1,iz,iq,3)=ddxN_q(it+1,iz,iq,3)+ddxN_z(it+1,iz,iq,3);
        ddxN(it+1,iz,iq,4)=ddxN_q(it+1,iz,iq,4)+ddxN_z(it+1,iz,iq,4);
        ddxN(it+1,iz,iq,5)=ddxN_q(it+1,iz,iq,5)+ddxN_z(it+1,iz,iq,5);
        
        U(it,iz,iq,1)=t_sample*uk_hat(it+1,iz,iq);
        U(it,iz,iq,2)=t_sample*uk_hat(it+2,iz,iq);
        U(it,iz,iq,3)=t_sample*uk_hat(it+3,iz,iq);
        U(it,iz,iq,4)=t_sample*uk_hat(it+4,iz,iq);
        U(it,iz,iq,5)=t_sample*uk_hat(it+5,iz,iq);
        U(it+1,iz,iq,1)=t_sample*uk_hat(it+2,iz,iq);
        U(it+1,iz,iq,2)=t_sample*uk_hat(it+3,iz,iq);
        U(it+1,iz,iq,3)=t_sample*uk_hat(it+4,iz,iq);
        U(it+1,iz,iq,4)=t_sample*uk_hat(it+5,iz,iq);
        U(it+1,iz,iq,5)=t_sample*uk_hat(it+6,iz,iq);
        
        dddxN_t(it,iz,iq,:)=(squeeze(ddxN(it+1,iz,iq,:))-squeeze(ddxN(it,iz,iq,:)))./t_sample;
        dU_t(it,iz,iq,:)=(squeeze(U(it+1,iz,iq,:))-squeeze(U(it,iz,iq,:)))./t_sample;      

        dxN_z(it,iz,iq,1)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+1,iz,iq));
        dxN_z(it,iz,iq,2)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+2,iz,iq));
        dxN_z(it,iz,iq,3)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+3,iz,iq));
        dxN_z(it,iz,iq,4)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+4,iz,iq));
        dxN_z(it,iz,iq,5)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+5,iz,iq));    
        dxN_q(it,iz,iq,1)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+1,iz,iq));
        dxN_q(it,iz,iq,2)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+2,iz,iq));
        dxN_q(it,iz,iq,3)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+3,iz,iq));
        dxN_q(it,iz,iq,4)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+4,iz,iq));
        dxN_q(it,iz,iq,5)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+5,iz,iq));
        
        dxN_z(it+1,iz,iq,1)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+2,iz,iq));
        dxN_z(it+1,iz,iq,2)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+3,iz,iq));
        dxN_z(it+1,iz,iq,3)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+4,iz,iq));
        dxN_z(it+1,iz,iq,4)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+5,iz,iq));
        dxN_z(it+1,iz,iq,5)=t_sample*F(1,1)*squeeze(dy1_hat_z(it+6,iz,iq));    
        dxN_q(it+1,iz,iq,1)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+2,iz,iq));
        dxN_q(it+1,iz,iq,2)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+3,iz,iq));
        dxN_q(it+1,iz,iq,3)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+4,iz,iq));
        dxN_q(it+1,iz,iq,4)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+5,iz,iq));
        dxN_q(it+1,iz,iq,5)=t_sample*F(1,2)*squeeze(dy1_hat_q(it+6,iz,iq));  
        
        dxN_zq(it,iz,iq,1)=dxN_z(it,iz,iq,1)+dxN_q(it,iz,iq,1);
        dxN_zq(it,iz,iq,2)=dxN_z(it,iz,iq,2)+dxN_q(it,iz,iq,2);
        dxN_zq(it,iz,iq,3)=dxN_z(it,iz,iq,3)+dxN_q(it,iz,iq,3);
        dxN_zq(it,iz,iq,4)=dxN_z(it,iz,iq,4)+dxN_q(it,iz,iq,4);
        dxN_zq(it,iz,iq,5)=dxN_z(it,iz,iq,5)+dxN_q(it,iz,iq,5);
        dxN_zq(it+1,iz,iq,1)=dxN_z(it+1,iz,iq,1)+dxN_q(it+1,iz,iq,1);
        dxN_zq(it+1,iz,iq,2)=dxN_z(it+1,iz,iq,2)+dxN_q(it+1,iz,iq,2);
        dxN_zq(it+1,iz,iq,3)=dxN_z(it+1,iz,iq,3)+dxN_q(it+1,iz,iq,3);
        dxN_zq(it+1,iz,iq,4)=dxN_z(it+1,iz,iq,4)+dxN_q(it+1,iz,iq,4);
        dxN_zq(it+1,iz,iq,5)=dxN_z(it+1,iz,iq,5)+dxN_q(it+1,iz,iq,5);

        ddxN_zq_t(it,iz,iq,:)=(squeeze(dxN_zq(it+1,iz,iq,:))-squeeze(dxN_zq(it,iz,iq,:)))./t_sample;
        end
      end

    for iz=2:N_z
        for iq=2:N_q
            for i=1:N
         dddxN_z(it,iz,iq,i)=(squeeze(ddxN(it,iz,iq,i))-squeeze(ddxN(it,iz-1,iq,i)))./z_sample;
         dddxN_q(it,iz,iq,i)=(squeeze(ddxN(it,iz,iq,i))-squeeze(ddxN(it,iz,iq-1,i)))./q_sample;
         dcF_z(it,iz,iq,i)=(squeeze(dxN_zq(it,iz,iq,i))-squeeze(dxN_zq(it,iz-1,iq,i)))./z_sample;
         dcF_q(it,iz,iq,i)=(squeeze(dxN_zq(it,iz,iq,i))-squeeze(dxN_zq(it,iz,iq-1,i)))./q_sample;
            end
        end
    end
     for iz=1:N_z-1
        for iq=1:N_q-1
            for i=1:N
         ddddxN_z(it,iz,iq,i)=(squeeze(dddxN_z(it,iz+1,iq,i))-squeeze(dddxN_z(it,iz,iq,i)))./z_sample;
         ddddxN_q(it,iz,iq,i)=(squeeze(dddxN_q(it,iz,iq+1,i))-squeeze(dddxN_q(it,iz,iq,i)))./q_sample;
         ddcF_z(it,iz,iq,i)=(squeeze(dcF_z(it,iz+1,iq,i))-squeeze(dcF_z(it,iz,iq,i)))./z_sample;
         ddcF_z(it,iz,iq,i)=(squeeze(dcF_z(it,iz,iq+1,i))-squeeze(dcF_z(it,iz,iq,i)))./q_sample;
            end
       end
     end 
     

        for iz=2:N_z-1  
           for iq=2:N_q-1 
        AF=AFlinear(Q*squeeze(xN(it,iz,iq,:))-Q*squeeze(dxN_zq(it,iz,iq,:))+Q*squeeze(ddxN(it,iz,iq,:))+Q*squeeze(U(it,iz,iq,:))); 
     
      xN(it+1,iz,iq,:)=squeeze(xN(it,iz,iq,:))+t_sample*Qpinv*(-Q*squeeze(-ddxN_zq_t(it,iz,iq,:))-Q*squeeze(dddxN_t(it,iz,iq,:))-Q*squeeze(dU_t(it,iz,iq,:))+D_hat*squeeze(DDxN(it,iz,iq,:))...
                       +D_hat*squeeze(ddddxN_z(it,iz,iq,:)+ddddxN_q(it,iz,iq,:))-D_hat*squeeze(ddcF_z(it,iz,iq,:)+ddcF_q(it,iz,iq,:))-gamma*AF);
   
          end
        end    

      A1(it+1)=2*a*a*(squeeze(xN(it+1,is1(1,1)-ia,is1(1,2)-ia,2))+squeeze(xN(it+1,is1(1,1)+ia,is1(1,2)-ia,2)))/2 ...
                +(squeeze(xN(it+1,is1(1,1)-ia,is1(1,2)+ia,2))+squeeze(xN(it+1,is1(1,1)+ia,is1(1,2)+ia,2)))/2;
      A2(it+1)=2*a*a*(squeeze(xN(it+1,is2(1,1)-ia,is2(1,2)-ia,2))+squeeze(xN(it+1,is2(1,1)+ia,is2(1,2)-ia,2)))/2 ...
                +(squeeze(xN(it+1,is2(1,1)-ia,is2(1,2)+ia,2))+squeeze(xN(it+1,is2(1,1)+ia,is2(1,2)+ia,2)))/2;
      A3(it+1)=2*a*a*(squeeze(xN(it+1,is3(1,1)-ia,is3(1,2)-ia,2))+squeeze(xN(it+1,is3(1,1)+ia,is3(1,2)-ia,2)))/2 ...
                +(squeeze(xN(it+1,is3(1,1)-ia,is3(1,2)+ia,2))+squeeze(xN(it+1,is3(1,1)+ia,is3(1,2)+ia,2)))/2;
      A4(it+1)=2*a*a*(squeeze(xN(it+1,is4(1,1)-ia,is4(1,2)-ia,2))+squeeze(xN(it+1,is4(1,1)+ia,is4(1,2)-ia,2)))/2 ...
                +(squeeze(xN(it+1,is4(1,1)-ia,is4(1,2)+ia,2))+squeeze(xN(it+1,is4(1,1)+ia,is4(1,2)+ia,2)))/2;
    A5(it+1)=2*a*a*(squeeze(xN(it+1,is5(1,1)-ia,is5(1,2)-ia,2))+squeeze(xN(it+1,is5(1,1)+ia,is5(1,2)-ia,2)))/2 ...
                +(squeeze(xN(it+1,is5(1,1)-ia,is5(1,2)+ia,2))+squeeze(xN(it+1,is5(1,1)+ia,is5(1,2)+ia,2)))/2;
ur1(it+1)=-k_u*A1(it+1);           
ur2(it+1)=-k_u*A2(it+1);
ur3(it+1)=-k_u*A3(it+1);
ur4(it+1)=-k_u*A4(it+1);   
ur5(it+1)=-k_u*A5(it+1);   

v1(it+1,1)=-ur1(it+1)*k_v*(squeeze(xN(it+1,is1(1,1)+ia,is1(1,2),2))-squeeze(xN(it+1,is1(1,1)-ia,is1(1,2),2)))/(2*a);
v1(it+1,2)=-ur1(it+1)*k_v*(squeeze(xN(it+1,is1(1,1),is1(1,2)+ia,2))-squeeze(xN(it+1,is1(1,1),is1(1,2)-ia,2)))/(2*a);

v2(it+1,1)=-ur2(it+1)*k_v*(squeeze(xN(it+1,is2(1,1)+ia,is2(1,2),2))-squeeze(xN(it+1,is2(1,1)-ia,is2(1,2),2)))/(2*a);
v2(it+1,2)=-ur2(it+1)*k_v*(squeeze(xN(it+1,is2(1,1),is2(1,2)+ia,2))-squeeze(xN(it+1,is2(1,1),is2(1,2)-ia,2)))/(2*a);

v3(it+1,1)=-ur3(it+1)*k_v*(squeeze(xN(it+1,is3(1,1)+ia,is3(1,2),2))-squeeze(xN(it+1,is3(1,1)-ia,is3(1,2),2)))/(2*a);
v3(it+1,2)=-ur3(it+1)*k_v*(squeeze(xN(it+1,is3(1,1),is3(1,2)+ia,2))-squeeze(xN(it+1,is3(1,1),is3(1,2)-ia,2)))/(2*a);

v4(it+1,1)=-ur4(it+1)*k_v*(squeeze(xN(it+1,is4(1,1)+ia,is4(1,2),2))-squeeze(xN(it+1,is4(1,1)-ia,is4(1,2),2)))/(2*a);
v4(it+1,2)=-ur4(it+1)*k_v*(squeeze(xN(it+1,is4(1,1),is4(1,2)+ia,2))-squeeze(xN(it+1,is4(1,1),is4(1,2)-ia,2)))/(2*a);
v5(it+1,1)=-ur5(it+1)*k_v*(squeeze(xN(it+1,is5(1,1)+ia,is5(1,2),2))-squeeze(xN(it+1,is5(1,1)-ia,is5(1,2),2)))/(2*a);
v5(it+1,2)=-ur5(it+1)*k_v*(squeeze(xN(it+1,is5(1,1),is5(1,2)+ia,2))-squeeze(xN(it+1,is5(1,1),is5(1,2)-ia,2)))/(2*a);

    s1(it+1,1)=s1(it,1)+t_sample*v1(it,1);
    s1(it+1,2)=s1(it,2)+t_sample*v1(it,2);
    s2(it+1,1)=s2(it,1)+t_sample*v2(it,1);
    s2(it+1,2)=s2(it,2)+t_sample*v2(it,2);
    s3(it+1,1)=s3(it,1)+t_sample*v3(it,1);
    s3(it+1,2)=s3(it,2)+t_sample*v3(it,2);
    s4(it+1,1)=s4(it,1)+t_sample*v4(it,1);
    s4(it+1,2)=s4(it,2)+t_sample*v4(it,2);
    s5(it+1,1)=s5(it,1)+t_sample*v5(it,1);
    s5(it+1,2)=s5(it,2)+t_sample*v5(it,2);
    
    is1=[round(s1(it+1,1)/z_sample) round(s1(it+1,2)/q_sample)];
    is2=[round(s2(it+1,1)/z_sample) round(s2(it+1,2)/q_sample)];
    is3=[round(s3(it+1,1)/z_sample) round(s3(it+1,2)/q_sample)];
    is4=[round(s4(it+1,1)/z_sample) round(s4(it+1,2)/q_sample)];
    is5=[round(s5(it+1,1)/z_sample) round(s5(it+1,2)/q_sample)];   
    
 for iz=1:N_z-1
  for iq=1:N_q-1         

    if is1(1,1)-ia<=iz && iz<=is1(1,1)+ia
        if is1(1,2)-ia<=iq && iq<=is1(1,2)+ia
        uk(it+1,iz,iq)=-k_u*A1(it+1);
        end
    end
    if is2(1,1)-ia<=iz && iz<=is2(1,1)+ia
        if is2(1,2)-ia<=iq && iq<=is2(1,2)+ia
        uk(it+1,iz,iq)=-k_u*A2(it+1);
        end
    end
    if is3(1,1)-ia<=iz && iz<=is3(1,1)+ia
        if is3(1,2)-ia<=iq && iq<=is3(1,2)+ia
        uk(it+1,iz,iq)=-k_u*A3(it+1);
        end
    end
    if is4(1,1)-ia<=iz && iz<=is4(1,1)+ia
        if is4(1,2)-ia<=iq && iq<=is4(1,2)+ia
        uk(it+1,iz,iq)=-k_u*A4(it+1);
        end
    end
    if is5(1,1)-ia<=iz && iz<=is5(1,1)+ia
        if is5(1,2)-ia<=iq && iq<=is5(1,2)+ia
        uk(it+1,iz,iq)=-k_u*A5(it+1);
        end
    end
    
  end     
 end
        
end

 y_y=zeros(N_z,N_q,N_t);
 put_uk=zeros(N_z,N_q,N_t);
 y1_y1=zeros(N_z,N_q,N_t);
 xN_xN=zeros(N_z,N_q,N_t);
 
 for it=1:N_t
 for iz=1:N_z
        for iq=1:N_q
            y_y(iz,iq,it)=y(it,iz,iq); 
            put_uk(iz,iq,it)=uk(it,iz,iq);
            y1_y1(iz,iq,it)=y1(it,iz,iq);
            xN_xN(iz,iq,it)=xN(it,iz,iq,2);
        end
 end
 end

Sum_y=zeros(N_t,N_z,N_q);
Sum_yy=zeros(N_t);
Sum_y1=zeros(N_t,N_z,N_q);
Sum_y1y1=zeros(N_t);
Sum_e=zeros(N_t,N_z,N_q);
Sum_ee=zeros(N_t);
Sum_xN1=zeros(N_t,N_z,N_q,N);
Sum_xNxN=zeros(N_t,N);
Sum_uk=zeros(N_t,N_z,N_q);
Sum_ukuk=zeros(N_t);
for it=1:N_t
    for iz=1:N_z-1
        for iq=1:N_q-1
          Sum_y(it,iz+1,iq+1)=Sum_y(it,iz,iq)+y(it,iz,iq)*y(it,iz,iq)*z_sample*q_sample;
          Sum_y1(it,iz+1,iq+1)=Sum_y1(it,iz,iq)+y1(it,iz,iq)*y1(it,iz,iq)*z_sample*q_sample;
          Sum_e(it,iz+1,iq+1)=Sum_e(it,iz,iq)+error(it,iz,iq)*error(it,iz,iq)*z_sample*q_sample;
          Sum_uk(it,iz+1,iq+1)=Sum_uk(it,iz,iq)+uk(it,iz,iq)*uk(it,iz,iq)*z_sample*q_sample;
          for i=1:N
          Sum_xN1(it,iz+1,iq+1,i)=Sum_xN1(it,iz,iq,i)+xN(it,iz,iq,i)*xN(it,iz,iq,i)*z_sample*q_sample;
          end
        end
    end
    Sum_yy(it)=Sum_y(it,N_z,N_q)^0.5;
%     Sum_yy_open(it)=Sum_y(it,N_z,N_q)^0.5;
    Sum_y1y1(it)=Sum_y1(it,N_z,N_q)^0.5;
    Sum_ee(it)=Sum_e(it,N_z,N_q)^0.5;
    Sum_ukuk(it)=Sum_uk(it,N_z,N_q)^0.5;
    for i=1:N
    Sum_xNxN(it,i)=Sum_xN1(it,N_z,N_q,i)^0.5;
    end
end

save RDNN_DSS_5R_2D_v1 ttt zzz qqq Sum_yy Sum_ee Sum_uk Sum_y1y1 Sum_xNxN put_uk y_y y y1 y1_hat N_t N_z N_q s1 s2 s3 s4 s5 xN uk u_o A1 A2 A3 A4 A5 A1_hat A2_hat A3_hat A4_hat A5_hat

% for it=1:N_t
%    if Sum_yy(it)<=0.05
%     T=(it+1)*0.01;
%         disp('The condition is satisfied and the loop is exited.');
%         fprintf('The convergence time is %f seconds。\n', T);
%         save PDNN_convergence_time T
%         return;
%    end
% end

% [min_value, min_index] = min(Sum_yy(:,1));

figure ()
plot(ttt,Sum_yy(:,1))
hold on
plot(ttt,Sum_y1y1(:,1))
grid on;
xlabel('t (sec.)');
ylabel('y(t)');
% xlim([0 5]);
h=legend('y','y1');
%%%%comparison
figure ()
plot(ttt,Sum_yy(:,1))
hold on
plot(ttt,Sum_yy_VPNN(:,1))
hold on
plot(ttt,Sum_yy_open(:,1))
grid on;
xlabel('t(sec.)');
ylabel('y(t)');
h=legend('PDNN','VPNN','open');

figure ()
plot(ttt,Sum_ee(:,1))
grid on;
xlabel('t');
ylabel('error');
% % xlim([0 5]);

figure()
for it=1:10:N_t
%     pause(0.2);
    fprintf('%d\n', it);
pcolor(zzz,qqq,y_y(:,:,it));hold on;
shading interp;
colorbar; colormap(jet);
plot(s1(it,2),s1(it,1),'*w','LineWidth', 3);hold on;
plot(s2(it,2),s2(it,1),'*m','LineWidth', 3);hold on;
plot(s3(it,2),s3(it,1),'*r','LineWidth', 3);hold on;
plot(s4(it,2),s4(it,1),'*c','LineWidth', 3);hold on;
plot(s5(it,2),s5(it,1),'*g','LineWidth', 3);
MakeGif('2D_RDNN_5robot.Gif',it);
xlabel('s2');ylabel('s1');
end

% for it=1:100:N_t
it=150;
figure()
set(gcf,'unit','centimeters','position',[10 5 13 11])
pcolor(zzz,qqq,y_y(:,:,it));hold on;
shading interp;
colorbar; colormap(jet);
plot(s1(it,2),s1(it,1),'*w','LineWidth', 5);hold on;
plot(s2(it,2),s2(it,1),'*m','LineWidth', 5);hold on;
plot(s3(it,2),s3(it,1),'*r','LineWidth', 5);hold on;
plot(s4(it,2),s4(it,1),'*c','LineWidth', 5);hold on;
plot(s5(it,2),s5(it,1),'*g','LineWidth', 5);
xlabel('s2');ylabel('s1');
caxis([0 1]);
title('t=1.5s');
% end

figure();
plot(s1(1,2),s1(1,1),'*k','LineWidth', 5);hold on;
plot(s2(1,2),s2(1,1),'*m','LineWidth', 5);hold on;
plot(s3(1,2),s3(1,1),'*r','LineWidth', 5);hold on;
plot(s4(1,2),s4(1,1),'*c','LineWidth', 5);hold on;
plot(s5(1,2),s5(1,1),'*g','LineWidth', 5);hold on;
plot(s1(:,2),s1(:,1),'k','LineWidth', 1);hold on;
plot(s2(:,2),s2(:,1),'m','LineWidth', 1);hold on;
plot(s3(:,2),s3(:,1),'r','LineWidth', 1);hold on;
plot(s4(:,2),s4(:,1),'c','LineWidth', 1);hold on;
plot(s5(:,2),s5(:,1),'g','LineWidth', 1);hold on;
xlabel('s2');
ylabel('s1');
xlim([0 3]);
ylim([0 3]);
title('Movement trajectories of robots');
h=legend('Robot1','Robot2','Robot3','Robot4','Robot5');

figure()
for it=1:10:N_t
    pause(0.1);
    fprintf('%d\n', it);
mesh(zzz,qqq,y_y(:,:,it));
xlabel('s2');
ylabel('s1');
zlabel('y(s1,s2)');
view(-37.5+90,20);
title('five robot control');
MakeGif('3D_RDNN_5R.Gif',it);
end

% it=1;
% figure()
% mesh(qqq,zzz,y_y(:,:,it));
% xlabel('$s_2$','interpreter','latex','FontSize',26);
% ylabel('$s_1$','interpreter','latex', 'FontSize',26);
% zlabel('$c(t,s_1,s_2)$','interpreter','latex','FontSize',26);
% zlim([0 1]);
% view(-37.5+90,20);
% title('$t$ = 0 s','interpreter','latex','FontSize',26);

