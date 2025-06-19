
load Data_writing_close
% load Data_writing_open

X1 = ZZ;                             
Y1 = QQ;                              
Z1 = YY;                       
Z2 = TT;                       

figure()
scatter3(X1,Y1,Z1,40,Z2,'filled')    
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('$s_2$','interpreter','latex','FontSize',18)
ylabel('$s_1$','interpreter','latex','FontSize',18)
zlabel('$c(t,s_1,s_2)$','interpreter','latex','FontSize',18);
colormap(jet);  

cb = colorbar; title(cb, '$t$ (s)','interpreter','latex','FontSize',18);
% cb.Label.String = 'Fatalities per 100M vehicle-miles';
