close all; clc; clear all;

% creaiamo una variabile struct per poi passarla come file json 
fusolage=struct;

% fattori di conversione
%in --> m 
in_m = 0.0254;
ft_m = 0.3048;
ft_in = 12;
lb_kg = 0.4536;

%passo parametri ala
load('..\wing DABCM98\caratteristiche_ala.mat')
b = b*ft_m;                     %m^2
S_wing = S_wing*ft_m^2;         %m^2




% parametri input
npax = 98;

seat_width = 20 * in_m;     %>16
fusolage.seat_width = seat_width;

seat_pitch = 30 * in_m;     %>30
fusolage.seat_picth = seat_pitch;

asile_width = 18 * in_m;    %>15
fusolage.asile_width = asile_width;

asile_height = 74 * in_m;   %>60
fusolage.asile_height = asile_height;

lavatories_pax = 2;         %1/50
fusolage.lavatories_pax  = lavatories_pax;

baggage_pax = 40 * lb_kg; 
fusolage.baggage_pax = baggage_pax; 

% uscite da normativa
typeI_width = 30 * in_m;    %da 24" a salire
ntypeI = 1;
fusolage.typeI_width = typeI_width;

typeIII_width = 20 * in_m;   %da 20"
ntypeIII = 2;                %per lato
fusolage.typeIII_width = typeIII_width;


% galley 1 
galley_width = 30 * in_m;   %da rivedere su tabelle ex dehavilland dhc-7
fusolage.galley_width = galley_width;


%toilets 1 per lato
toilet_width = 36 * in_m;   %
fusolage.toilet_width = toilet_width;


%scelta confiurazione full_economy con  
% 1)2+2 
% 2)3+2
conf = 2;

%dimensioni sezione fusoliera e numero di file (approssimano per eccesso)
switch conf
    case  1
       
        cabin_width = 4 * seat_width + asile_width;
        % parte strutturale 7% della larghezza
        struttura =0.07 * cabin_width;
        cabin_width = cabin_width + 0.08 * cabin_width;
        nfile = ceil(npax/4);
        
        fusolage.cabin_width = cabin_width;
        fusolage.struttura  = struttura;
        fusolage.nfile = nfile;
      
        
        
    case 2
        
        cabin_width = 5 * seat_width + asile_width;
        % parte strutturale 7% della larghezza
        struttura =0.07 * cabin_width;
        cabin_width = cabin_width + 0.07 * cabin_width;
        nfile = ceil(npax/5);
          
        fusolage.cabin_width = cabin_width;
        fusolage.struttura  = struttura;
        fusolage.nfile = nfile;
        
end

df = cabin_width;


% rapporto finezza nose da 1.2 a 2.5
ln_df = 1.2;   
fusolage.ln_df = ln_df;

ln = ln_df * df;
fusolage.ln = ln;


% rapporto finezza tailcone da 2 a 5 , tipico 2.8-3.2
lt_df = 2.8;
fusolage.lf_df = lt_df;

lt = lt_df * df;
fusolage.lt = lt;

% lunghezza corpo fusoliera(cabina) (non considero la galley perche ho approssimato
% per eccesso il numero di file)
lc = nfile*seat_pitch + ntypeIII*typeIII_width + ntypeI*typeI_width + toilet_width;
fusolage.lc = lc;

% lunghezza complessiva della fusoliera

lf = lc + lt +ln;
fusolage.lf = lf;

% rapporto di finezza 
lf_df = lf/df;
fusolage.lf_df = lf_df;


%angolo upsweep 
teta_f = 16;    %deg
fusolage.teta_f = teta_f;

%Windshield angle
psi = 45;       %deg
fusolage.psi = psi;


%% fusdes method
%dalla geometria definita
% da aerei simili l'ala sarà posizionata circa al 0.43 di lf
% per cui xref=0.5 --->posizione ripetto cg
S_front = pi*df^2/4;        %m^2

S_wet_nose = pi*df*0.75*ln; %m^2
S_wet_cabin = pi*df*lc;     %m^2
S_wet_tail = pi*df*0.72*lt; %m^2

S_wet = S_wet_nose + S_wet_cabin + S_wet_tail; %m^2

Kn = 1.95;
Kc = 1.08;
Kt = 0.73;

% in cruise
Re_fus = 0.549*170*lf/1.55e-5;
M_cr = 0.55;

% coeff resistenza
C_D_fp = 0.455/(log(Re_fus)^2.58+(1+0.144*M_cr^2))^0.65;
C_D_fus__ = (Kn*S_wet_nose/S_wet + Kc*S_wet_cabin/S_wet +...
    Kt*S_wet_tail/S_wet)*C_D_fp*S_wet/S_front;


CD_fus = C_D_fus__*S_front/S_wing


% coefficiente di momento
MAC = MAC * ft_m;
CM0fr = -0.30 - (-0.30 + 0.26)*7/15;   % xref =0.5*lf
D_CM0nose = -0.08; 
D_CM0tail = 0.0125;

CM0_fus__ = CM0fr + D_CM0nose + D_CM0tail; 
CM0_fus = CM0_fus__*S_front/S_wing*df/MAC


CM_a_fr = 0.22 - (0.22-0.195)*7/15; % xref =0.5*lf
CM_a_nose = -0.0035;
CM_a_tail = -0.002;
CM_a_fus__ =  CM_a_fr + CM_a_nose + CM_a_tail;

CM_a_fus = CM_a_fus__*S_front/S_wing*df/MAC

% plot del contributo legato alla fusoliera
% alpha = linspace(-5,10,16);
% Cm = CM0_fus + CM_a_fus*alpha;
% plot(alpha,Cm,'k')

% coefficiente  di yaw


CN_b_fr = -0.212 - (-0.212 + 0.190)*7/15; % xref =0.42*lf
CN_b_nose = +0.0023 - (0.0023-0.0032)*7/15;
CN_b_tail = 0-(0.-0.0008)*7/15;
CN_b_fus__ =  CN_b_fr + CN_b_nose + CN_b_tail;
CN_b_fus_isolated = CN_b_fus__*S_front/S_wing*df/b;



% ottenere il dato finale stabiliamo quella che sarà la configurazione
% del nostro aereo 
% configurazione high wing ---> z
Kwf = 1.16;
Kv = 0.87;
Kh = 0.975;

CN_b_fus = Kwf*Kv*Kh*CN_b_fus_isolated
%%

% json = jsonencode(fusolage)
% fid=fopen('fusolage.json','w');
% fprintf(fid,json);
% fclose('all')