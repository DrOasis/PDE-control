%%%Open loop
%%%2D space£¬5 robots

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


t0=0;
tf=10;
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
y1=zeros(N_t,N_z,N_q);
dy1_z=zeros(N_t,N_z,N_q);
ddy1_z=zeros(N_t,N_z,N_q);
dy1_q=zeros(N_t,N_z,N_q);
ddy1_q=zeros(N_t,N_z,N_q);

yout=zeros(N_t,N_z,N_q);
y1out=zeros(N_t,N_z,N_q);
u_o=zeros(N_t,N_z,N_q);

error=zeros(N_t,N_z,N_q);

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
% y1(1,iz,iq)=rand;
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

for it=1:N_t-1
    it  
    for iz=2:N_z-1  
        for iq=2:N_q-1
         y(it+1,iz,iq)=squeeze(y(it,iz,iq))+t_sample*(squeeze(D(1,1)*ddy_z(it,iz,iq)+D(2,2)*ddy_q(it,iz,iq))-squeeze(F(1,1)*dy_z(it,iz,iq)+F(1,2)*dy_q(it,iz,iq)));

         y1(it+1,iz,iq)=squeeze(y1(it,iz,iq))+t_sample*(squeeze(D(1,1)*ddy1_z(it,iz,iq)+D(2,2)*ddy1_q(it,iz,iq))-squeeze(F(1,1)*dy1_z(it,iz,iq)...
                      +F(1,2)*dy1_q(it,iz,iq))+squeeze(u_o(it,iz,iq)));
                  
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
         end
     end
    
end

% for it=1:N_t
%     for iz=1:N_z
%         y_l1(it,iz)=y(it,iz,l1);
%         output_u_o1(it,iz)=u_o(it,iz,l1);
%     end
% end

 y_y=zeros(N_z,N_q,N_t);
 
 for it=1:N_t
 for iz=1:N_z
        for iq=1:N_q
            y_y(iz,iq,it)=y(it,iz,iq);
        end
 end
 end

Sum_y=zeros(N_t,N_z,N_q);
Sum_yy_open=zeros(N_t);
Sum_y1=zeros(N_t,N_z,N_q);
Sum_y1y1=zeros(N_t);
Sum_e=zeros(N_t,N_z,N_q);
Sum_ee=zeros(N_t);
% Sum_y2=zeros(N_t,N_z,N_q);
for it=1:N_t
    for iz=1:N_z-1
        for iq=1:N_q-1
          Sum_y(it,iz+1,iq+1)=Sum_y(it,iz,iq)+y(it,iz,iq)*y(it,iz,iq)*z_sample*q_sample;
          Sum_y1(it,iz+1,iq+1)=Sum_y1(it,iz,iq)+y1(it,iz,iq)*y1(it,iz,iq)*z_sample*q_sample;
%  
%           Sum_y2(it,iz+1,iq+1)=Sum_y2(it,iz,iq)+y2(it,iz,iq)*y2(it,iz,iq)*z_sample*q_sample;
          
          Sum_e(it,iz+1,iq+1)=Sum_e(it,iz,iq)+error(it,iz,iq)*error(it,iz,iq)*z_sample*q_sample;

        end
    end
    Sum_yy_open(it)=Sum_y(it,N_z,N_q)^0.5;
    Sum_y1y1(it)=Sum_y1(it,N_z,N_q)^0.5;
    Sum_ee(it)=Sum_e(it,N_z,N_q)^0.5;
%     Sum_y2y2(it)=Sum_y2(it,N_z,N_q)^0.5;
end

% save Open_loop y_y y1 N_t ttt zzz qqq Sum_yy_open Sum_ee Sum_y1y1 

figure ()
plot(ttt,Sum_yy_open(:,1))
hold on
plot(ttt,Sum_y1y1(:,1))
grid on;
xlabel('t (sec.)');
ylabel('y(t)');
% xlim([0 5])
h=legend('y','y1');

figure ()
plot(ttt,Sum_ee(:,1))
grid on;
xlabel('t (sec.)');
ylabel('e(t)');



