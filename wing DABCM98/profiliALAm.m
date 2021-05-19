close all; clear all; clc;
N618=readmatrix('dati profili\NACA63618.dat');
N615=readmatrix('dati profili\NACA63615.dat');


figure
plot(N618(:,1),N618(:,2),'k'); hold on;
axis image; axis([-0.05 1.05 -0.2 0.2]);
xlabel('$\frac{x}{c}$','Interpreter','latex','FontName','Palatino')
ylabel('$\frac{z}{c}$','Interpreter','latex','FontName','Palatino');
    
legend({'Profilo 63-618'},'Location','northeast','Interpreter','latex','FontName','Palatino')
saveas(gcf,'risultati immagini\profilo_63618.png')

figure
plot(N615(:,1),N615(:,2),'k'); hold on;
axis image; axis([-0.05 1.05 -0.2 0.2]);
xlabel('$\frac{x}{c}$','Interpreter','latex','FontName','Palatino')
ylabel('$\frac{z}{c}$','Interpreter','latex','FontName','Palatino');
    
legend({'Profilo 63-615'},'Location','northeast','Interpreter','latex','FontName','Palatino')
saveas(gcf,'risultati immagini\profilo_63615.png')

figure
plot(N618(:,1),N618(:,2),'k',N615(:,1),N615(:,2),'-.k');  
axis image; 
axis([-0.05 1.05 -0.2 0.2]);
xlabel('$\frac{x}{c}$','Interpreter','latex','FontName','Palatino')
ylabel('$\frac{z}{c}$','Interpreter','latex','FontName','Palatino');
legend({'Profilo 63-618','Profilo 63-615'},'Location','northeast','Interpreter','latex','FontName','Palatino')
saveas(gcf,'risultati immagini\confronto_profili_618_615.png')

%%
Q1 = readmatrix('polari profili\polari_NACA63615_Re6p5e6.dat');
Q2 = readmatrix('polari profili\polari_NACA63615_Re11e6.dat');
R1 = readmatrix('polari profili\polari_NACA63618_Re12p5e6.dat');
R2 = readmatrix('polari profili\polari_NACA63618_Re19e6.dat');


NACA63615Re6p5e6  = Q1(1:end,:);
NACA63615Re11e6 = Q2(1:end,:);
NACA63618Re12p5e6  = R1(1:end,:);
NACA63618Re19e6 = R2(1:end,:);


a1=NACA63615Re6p5e6(:,1);
a2=NACA63615Re11e6(:,1);
a3=NACA63618Re12p5e6(:,1);
a4=NACA63618Re19e6(:,1);

Cl1=NACA63615Re6p5e6(:,2);
Cl2=NACA63615Re11e6(:,2);
Cl3=NACA63618Re12p5e6(:,2);
Cl4=NACA63618Re19e6(:,2);


Cd1=NACA63615Re6p5e6(:,3);
Cd2=NACA63615Re11e6(:,3);
Cd3=NACA63618Re12p5e6(:,3);
Cd4=NACA63618Re19e6(:,3);


Cm1=NACA63615Re6p5e6(:,4);
Cm2=NACA63615Re11e6(:,4);
Cm3=NACA63618Re12p5e6(:,4);
Cm4=NACA63618Re19e6(:,4);


figure 
plot(a1,Cl1,'-k');
xline(0,'k');
yline(0,'k');
axis([-12 35 -1 2.4]);
legend({'NACA 63615 $R_e=6.50\times10^{6}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
xlabel('\textbf{$\alpha^\circ$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
grid on;
cl1max=max(Cl1)
saveas(gcf,'risultati immagini\curvacl_615_Re6p5e6.png')


figure 
plot(a2,Cl2,'-k');
xline(0,'k');
yline(0,'k');
axis([-12 35 -1 2.4]);
legend({'NACA 63615 $R_e=1.10\times10^{7}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
xlabel('\textbf{$\alpha^\circ$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
grid on;
cl2max=max(Cl2)
saveas(gcf,'risultati immagini\curvacl_615_Re11e6.png')


figure 
plot(a3,Cl3,'-k');
xline(0,'k');
yline(0,'k');
axis([-12 35 -1 2.4]);
legend({'NACA 63618 $R_e=1.25\times10^{7}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
xlabel('\textbf{$\alpha^\circ$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
grid on;
cl3max=max(Cl3)
saveas(gcf,'risultati immagini\curvacl_18_Re12p5e6.png')


figure 
plot(a4,Cl4,'-k');
xline(0,'k');
yline(0,'k');
axis([-12 35 -1 2.4]);
legend({'NACA 63618 $R_e=1.90\times10^{7}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
xlabel('\textbf{$\alpha^\circ$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
grid on;
cl4max=max(Cl4)
saveas(gcf,'risultati immagini\curvacl_618_Re19e6.png')

%%
P1 = readmatrix('polari profili\polari_NACA63615_eulero.dat');
P2 = readmatrix('polari profili\polari_NACA63618_eulero.dat');
Cm15 = P1(:,4);
Cm18 = P2(:,4);
Cl15 = P1(:,2);
Cl18 = P2(:,2);
alpha15 = P1(:,1);
alpha18 = P2(:,1);


xac15=.25-(Cm15(4)-Cm15(5))/Cl15(4)
xac18=.25-(Cm18(4)-Cm18(5))/Cl18(4)

Clalpha15EU=(Cl15(5)-Cl15(4))/(alpha15(5)-alpha15(4))
Clalpha18EU=(Cl18(5)-Cl18(4))/(alpha18(5)-alpha18(4))

%%
figure 
plot(Cd1*10000,Cl1,'k',Cd2*10000,Cl2,'-.k'); 
axis([0 1000 -1 2.4]);
lgd = legend({'$R_e=6.50\times10^{6}$','$R_e=1.10\times10^{7}$'},'Location','east','Interpreter','latex','FontName','Palatino');
title(lgd,'NACA 63-615','Interpreter','latex','FontName','Palatino')
xlabel('\textbf{$C_d$ drag count}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
%text(800,-0.5,'NACA 63615','Interpreter','latex','FontName','Palatino')
saveas(gcf,'risultati immagini\polare_resistenza615.png')


figure 
plot(Cd3*10000,Cl3,'k',Cd2*10000,Cl2,'-.k'); 
axis([0 1000 -1 2.4]);
lgd = legend({'$R_e=1.25\times10^{7}$','$R_e=1.90\times10^{7}$'},'Location','east','Interpreter','latex','FontName','Palatino');
title(lgd,'NACA 63-618','Interpreter','latex','FontName','Palatino')
xlabel('\textbf{$C_d$ drag count}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
%text(800,-0.5,'NACA 63618','Interpreter','latex','FontName','Palatino')
saveas(gcf,'risultati immagini\polare_resistenza618.png')

%% coefficienti di momento

figure 
plot(Cl2,Cm2,'k'); 
axis([-1 2.4 -0.35 0.05]);
xlabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_{m_{\frac{c}{4}}}$}','Interpreter','latex','FontName','Palatino');xline(0,'k');
yline(0,'k');
legend({'NACA 63615 $R_e=1.10\times10^{6}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
yticks(-.35:0.05:0.05);
pbaspect([3 1 1])
saveas(gcf,'risultati immagini\cmc4_615_Re11e6.png')


figure 
plot(Cl4,Cm4,'k'); 
axis([-1 2.4 -0.35 0.05]);
xlabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_{m_{\frac{c}{4}}}$}','Interpreter','latex','FontName','Palatino');xline(0,'k');
yline(0,'k');
legend({'NACA 63618 $R_e=1.90\times10^{7}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
yticks(-.35:0.05:-0.05);
pbaspect([3 1 1])
saveas(gcf,'risultati immagini\cmc4_618_Re19e6.png')


%Coefficiente di momento rispetto al centro aerodinamico
for i=1:length(Cl2)
Cmac2(i)= Cm2(i)+Cl2(i)*(xac15-.25);
end
figure 
plot(Cl2,Cmac2,'k'); 
axis([-1 2.4 -0.35 0.05]);
xlabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_{m,ac}$}','Interpreter','latex','FontName','Palatino');xline(0,'k');
yline(0,'k');
legend({'NACA 63615 $R_e=1.10\times10^{6}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
yticks(-.35:0.05:0.05);
pbaspect([3 1 1])
saveas(gcf,'risultati immagini\cmac_615_Re11e6.png')



for i=1:length(Cl4)
Cmac4(i)= Cm4(i)+Cl4(i)*(xac18-.25);
end
figure 
plot(Cl4,Cmac4,'k'); 
axis([-1 2.4 -0.35 0.05]);
xlabel('\textbf{$C_l$}','Interpreter','latex','FontName','Palatino');
ylabel('\textbf{$C_{m,ac}$}','Interpreter','latex','FontName','Palatino');xline(0,'k');
yline(0,'k');
legend({'NACA 63618 $R_e=1.90\times10^{7}$'},'Location','northwest','Interpreter','latex','FontName','Palatino');
yticks(-.35:0.05:-0.05);
pbaspect([3 1 1])
saveas(gcf,'risultati immagini\cmac_618_Re19e6.png')

%% 
A1=NACA63615Re6p5e6(:,1:4);
A2=NACA63615Re11e6(:,1:4);
A3=NACA63618Re12p5e6(:,1:4);
A4=NACA63618Re19e6(:,1:4);

save('polari profili\Naca615_0fl.dat','A1','-ascii');
save('polari profili\Naca615_250fl.dat','A2','-ascii');
save('polari profili\Naca618_0fl.dat', 'A3','-ascii');
save('polari profili\Naca618_250fl.dat', 'A4','-ascii');

%% 2d output
prof2d=struct
prof2d.Cla18 = Clalpha18;
prof2d.Cla15 = Clalpha15;
prof2d.Xac15 = xac15;
prof2d.Xac18 = xac18;
prof2d.Clmax13 = cl2max;
prof2d.Clmax18 = cl4max;

json = jsonencode(prof2d)
fid=fopen('profili.json','w');
fprintf(fid,json);
fclose('all')