close all; clc; clear all;
% creaiamo una variabile struct per poi passarla come file json 
fusolage=struct;

% fattori di conversione
%in --> m 
in_m = 0.0254;
ft_m = 0.3048;
ft_in = 12;
lb_kg = 0.4536;

 
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
ln_df = 1.2   
fusolage.ln_df = ln_df;

ln = ln_df * df;
fusolage.ln = ln;


% rapporto finezza tailcone da 2 a 5 , tipico 2.8-3.2
lt_df = 2.8;
fusolage.lf_df = lt_df;

lt = lt_df * df;
fusolage.lt = lt;

% lunghezza corpo fusoliera (non considero la galley perche ho approssimato
% per eccesso il numero di file
lc = nfile*seat_pitch + ntypeIII*typeIII_width + ntypeI*typeI_width + toilet_width;
fusolage.lc = lc;

% lunghezza complessiva della fusoliera

lf = lc + lt +ln;
fusolage.lf = lf;

% rapporto di finezza 
lf_df = lf/df;
fusolage.lf_df = lf_df;


%angolo upsweep 
teta_f = 20; %deg
fusolage.teta_f = teta_f;

json = jsonencode(fusolage)
fid=fopen('fusolage.json','w');
fprintf(fid,json);
fclose('all')


