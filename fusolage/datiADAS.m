close all; clc; clear all;

a=readmatrix('D:\AIRCRAFT DESIGN\From github\Aircraft-Design\fusolage\from ADAS\polariFusoliera.dat');

% CD
figure 
plot(a(:,1),a(:,8),'k');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$CD_{fus}$','Interpreter','latex');
grid minor;
% saveas(gcf,'D:\AIRCRAFT DESIGN\From github\Aircraft-Design\fusolage\immaginiFusoliera\CD_fusoliera.png');

% CM
figure 
plot(a(:,1),a(:,9),'k');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$CM_{fus}$','Interpreter','latex');
text(0,0.3,'','Interpreter','latex');
xline(0)
yline(0)
grid minor
%saveas(gcf,'D:\AIRCRAFT DESIGN\From github\Aircraft-Design\fusolage\immaginiFusoliera\CM_fusoliera.png');
% matlab2tikz('D:\AIRCRAFT DESIGN\From github\Aircraft-Design\fusolage\immaginiFusoliera\CM_fusoliera.tex')

