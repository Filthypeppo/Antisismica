clear all
clc
%%%NEWMARK%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULO PARA LA COMPONENTE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//90º//%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


 T =0.02:0.01:10; %vector de periodos (seg)
 maxu = zeros(length(T),1);
 maxv = zeros(length(T),1);
 maxa= zeros(length(T),1);
 
 
 m = 1; %ton
 load Comp90.txt 
 %%aguante = Sismo(255:length(Sismo),1);
 t = Comp90(:,1); % importo el vector de tiempo [s]
 Sv = zeros(length(T),3); %% creo el vector de pseudovelocidades.
 as = Comp90(:,2); % importo el vector de aceleraciones [g]
 as90 = as; %% usado para grafico linea base
 asglobal(:,1) = as; 
 p = -m*as*981; % tonf*cm/s2 ( vector de fuerzas ) 
 pt = zeros(length(t),1);
 dt = t(2); %extraigo la segunda fila del vector de tiempos el cual corresponde al Delta t
 var =1;
for Te = T  %recorremos el vector T
   
    w = 2*pi/Te; %rad/s
    k = (w^2)*m; % ton/(s^2)
    xi =0.05;
    c = 2*m*w*xi; % ton/s
    
    %w = sqrt(k*981/m); %rad/s

    %plot(t,as)

    u = zeros(length(t),1); %  u(1) = u0
    a = zeros(length(t),1);  %  a(1) = a0
    v = zeros(length(t),1);  %  v(1) = v0

    beta = 1/6;
    gamma = 1/2;

    a1 = ((1/(beta*(dt^2)))*m) + ((gamma/(beta*dt))*c); %ton/s2
    a2 = (1/(beta*(dt))*m) + ((gamma/beta)-1)*c; %ton/s
    a3 = (((1/(2*beta))-1)*m) + (dt*c*((gamma/(2*beta))-1));

    kt = k+a1; %ton/s^2

    for i=2:length(u)

        pt(i,1) = p(i,1) + a1*u(i-1,1) + a2*v(i-1,1)+ a3*a(i-1,1); %ton *cm/s2
        u(i,1) = pt(i,1)/kt; %cm
        v(i,1) = ((gamma/(beta*dt))*(u(i,1)-u(i-1,1)))+((1-(gamma/beta))*v(i-1,1)) + (dt*a(i-1,1)*(1-(gamma/(2*beta)))); %cm/s
        a(i,1) = ((1/(beta*(dt^2)))*(u(i,1)-u(i-1,1))) - (1/(beta*dt)*v(i-1,1)) - (((1/(2*beta))-1)*a(i-1,1)); %cm/s2;
    end

    maxu(var) =max(abs(u)); %cm 
    maxv(var) = max(abs(v)); %cm/s
    maxa(var) = max(abs(a))/981;  %(cm/s2)/g
    
    var = var+1;
end
Tplot = transpose(T);
w = 2*pi./Tplot;
Sv(:,1) = maxu.*w; % Pseudo velocidades para la componente de 90º
Sa = (maxu.*(w.^2))/981;
disp(['PGA Componente 90 app = ',num2str(Sa(1))])
%subplot(1,2,1),plot(T,maxa);
subplot(1,3,1),plot(T,Sa,'b');
xlim([0 4])
xlabel('Periodo [s]');
ylabel('Pseudoaceleraciones [g]');
title(['Componente 90º  PGA = ' , num2str(Sa(1)),' [g]'])
hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULO PARA LA COMPONENTE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//360º//%%%%%%%%%%%%%%%%%%%%%%%%%%



 

%%%NEWMARK%%%

 T =0.02:0.01:10; %vector de periodos (seg)
 maxu = zeros(length(T),1);
 maxv = zeros(length(T),1);
 maxa= zeros(length(T),1);
 
 m = 1; %ton
 load Comp360.txt 
 %%aguante = Sismo(255:length(Sismo),1);
 t = Comp360(:,1); % importo el vector de tiempo [s]
 as = Comp360(:,2); % importo el vector de aceleraciones [g]
 as360 = as; %% usado para grafico linea base
 p = -m*as*981; % tonf*cm/s2 ( vector de fuerzas ) 
 pt = zeros(length(t),1);
 dt = t(2); %extraigo la segunda fila del vector de tiempos el cual corresponde al Delta t
 var =1;
for Te = T  %recorremos el vector T
   
    w = 2*pi/Te; %rad/s
    k = (w^2)*m; % ton/(s^2)
    xi =0.05;
    c = 2*m*w*xi; % ton/s
    
    %w = sqrt(k*981/m); %rad/s

    %plot(t,as)

    u = zeros(length(t),1); %  u(1) = u0
    a = zeros(length(t),1);  %  a(1) = a0
    v = zeros(length(t),1);  %  v(1) = v0

    beta = 1/6;
    gamma = 1/2;

    a1 = ((1/(beta*(dt^2)))*m) + ((gamma/(beta*dt))*c); %ton/s2
    a2 = (1/(beta*(dt))*m) + ((gamma/beta)-1)*c; %ton/s
    a3 = (((1/(2*beta))-1)*m) + (dt*c*((gamma/(2*beta))-1));

    kt = k+a1; %ton/s^2

    for i=2:length(u)

        pt(i,1) = p(i,1) + a1*u(i-1,1) + a2*v(i-1,1)+ a3*a(i-1,1); %ton *cm/s2
        u(i,1) = pt(i,1)/kt; %cm
        v(i,1) = ((gamma/(beta*dt))*(u(i,1)-u(i-1,1)))+((1-(gamma/beta))*v(i-1,1)) + (dt*a(i-1,1)*(1-(gamma/(2*beta)))); %cm/s
        a(i,1) = ((1/(beta*(dt^2)))*(u(i,1)-u(i-1,1))) - (1/(beta*dt)*v(i-1,1)) - (((1/(2*beta))-1)*a(i-1,1)); %cm/s2;
    end

    maxu(var) =max(abs(u)); %cm 
    maxv(var) = max(abs(v)); %cm/s
    maxa(var) = max(abs(a))/981;  %(cm/s2)/g
    
    var = var+1;
end
Tplot = transpose(T);
w = 2*pi./Tplot;
Sv(:,2) = maxu.*w; %% Pseudo velocidades para la componente 360º
Sa = (maxu.*(w.^2))/981;
disp(['PGA Componente 360 app = ',num2str(Sa(1))])
%subplot(1,2,1),plot(T,maxa);
subplot(1,3,2),plot(T,Sa,'r');
xlim([0 4])
xlabel('Periodo [s]');
ylabel('Pseudoaceleraciones [g]');
title(['Componente 360º  PGA = ' , num2str(Sa(1)),' [g]'])
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULO PARA LA COMPONENTE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//upº//%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%NEWMARK%%%

 T =0.02:0.01:10; %vector de periodos (seg)
 maxu = zeros(length(T),1);
 maxv = zeros(length(T),1);
 maxa= zeros(length(T),1);
 
 m = 1; %ton
 load Compup.txt 
 %%aguante = Sismo(255:length(Sismo),1);
 t = Compup(:,1); % importo el vector de tiempo [s]
 as = Compup(:,2); % importo el vector de aceleraciones [g]
 asup = as; %% usado para grafico linea base
 p = -m*as*981; % tonf*cm/s2 ( vector de fuerzas ) 
 pt = zeros(length(T),1);
 dt = t(2); %extraigo la segunda fila del vector de tiempos el cual corresponde al Delta t
 var =1;
for Te = T  %recorremos el vector T
   
    w = 2*pi/Te; %rad/s
    k = (w^2)*m; % ton/(s^2)
    xi =0.05;
    c = 2*m*w*xi; % ton/s
    
    %w = sqrt(k*981/m); %rad/s

    %plot(t,as)

    u = zeros(length(t),1); %  u(1) = u0
    a = zeros(length(t),1);  %  a(1) = a0
    v = zeros(length(t),1);  %  v(1) = v0

    beta = 1/6;
    gamma = 1/2;

    a1 = ((1/(beta*(dt^2)))*m) + ((gamma/(beta*dt))*c); %ton/s2
    a2 = (1/(beta*(dt))*m) + ((gamma/beta)-1)*c; %ton/s
    a3 = (((1/(2*beta))-1)*m) + (dt*c*((gamma/(2*beta))-1));

    kt = k+a1; %ton/s^2

    for i=2:length(u)

        pt(i,1) = p(i,1) + a1*u(i-1,1) + a2*v(i-1,1)+ a3*a(i-1,1); %ton *cm/s2
        u(i,1) = pt(i,1)/kt; %cm
        v(i,1) = ((gamma/(beta*dt))*(u(i,1)-u(i-1,1)))+((1-(gamma/beta))*v(i-1,1)) + (dt*a(i-1,1)*(1-(gamma/(2*beta)))); %cm/s
        a(i,1) = ((1/(beta*(dt^2)))*(u(i,1)-u(i-1,1))) - (1/(beta*dt)*v(i-1,1)) - (((1/(2*beta))-1)*a(i-1,1)); %cm/s2;
    end

    maxu(var) =max(abs(u)); %cm 
    maxv(var) = max(abs(v)); %cm/s
    maxa(var) = max(abs(a))/981;  %(cm/s2)/g
    
    var = var+1;
end
Tplot = transpose(T);
w = 2*pi./Tplot;
Sv(:,3) = maxu.*w;
Sa = (maxu.*(w.^2))/981;
%subplot(1,2,1),plot(T,maxa);
disp(['PGA Componente up app = ',num2str(Sa(1))])
subplot(1,3,3),plot(T,Sa,'g');
xlim([0 4])
xlabel('Periodo [s]');
ylabel('Pseudoaceleraciones [g]');
title(['Componente up  PGA = ' , num2str(Sa(1)),' [g]'])



% Obtengo la PGA, PGV y PGD en fundcion de las pseudo aceleraciones
% utilizando un metodo de integracion.

PGA90 = max(as90)*981;
PGA360 = max(as360)*981;
PGAUP = max(asup)*981;

vs90 = as90*0;
vs360 = as360*0;
vsup = asup*0;

%interal de aceleraciones ( VELOCIDADES )

for i =2:length(vs90)
    rec90 = dt*((as90(i) + as90(i-1))/2)*981; 
    vs90(i) = rec90 + vs90(i-1);
end
for i =2:length(vs360)
    rec360 = dt*((as360(i) + as360(i-1))/2)*981; 
    vs360(i) = rec360 + vs360(i-1);
end
for i =2:length(vsup)
    recup = dt*((asup(i) + asup(i-1))/2)*981; 
    vsup(i) = recup + vsup(i-1);
end

PGV90 = max(vs90);
PGV360 = max(vs360);
PGVUP = max(vsup);


%%%%%% INTEGRAL DE VELOCIDADES %%%% ( DESPLAZAMIENTOS)

us90 = vs90*0;
us360 = vs360*0;
usup = vsup*0;


for i =2:length(us90)
    rec90 = dt*((vs90(i) + vs90(i-1))/2); 
    us90(i) = rec90 + us90(i-1);
end
for i =2:length(us360)
    rec360 = dt*((vs360(i) + vs360(i-1))/2); 
    us360(i) = rec360 + us360(i-1);
end
for i =2:length(usup)
    recup = dt*((vsup(i) + vsup(i-1))/2); 
    usup(i) = recup + usup(i-1);
end


PGD90 = max(us90);
PGD360 = max(us360);
PGDUP = max(usup);

Ttetra =0.02:0.01:10;

PGA90vec = ones(1,length(Ttetra))*PGA90./w; %cm/s
PGA360vec = ones(1,length(Ttetra))*PGA360./w; %cm/s
PGAUPvec = ones(1,length(Ttetra))*PGAUP./w; %cm/s

PGV90vec = ones(1,length(Ttetra))*PGV90; %cm/s
PGV360vec = ones(1,length(Ttetra))*PGV360; %cm/s
PGVUPvec = ones(1,length(Ttetra))*PGVUP; %cm/s


PGD90vec = ones(1,length(Ttetra))*PGD90.*w; %cm/s
PGD360vec = ones(1,length(Ttetra))*PGD360.*w; %cm/s
PGDUPvec = ones(1,length(Ttetra))*PGDUP.*w; %cm/s

LB90 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 
LB360 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 
LBUP = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 

for i =1:length(Ttetra)
   LB90(i) = min([PGA90vec(i) PGV90vec(i) PGD90vec(i)]);
end

for i =1:length(Ttetra)
   LB360(i) = min([PGA360vec(i) PGV360vec(i) PGD360vec(i)]);
end

for i =1:length(Ttetra)
   LBUP(i) = min([PGAUPvec(i) PGVUPvec(i) PGDUPvec(i)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Factores de amplificacion de espectros de diseño elastico para
%%xi = 5%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alfa50 = ([2.12 1.65 1.39]);
alfa84 = ([2.71 2.30 2.01]);

%AMPLIFICACION PARA CADA COMPONENTE:

%COMPONENTE 90º
PGA90vec50 = PGA90vec*alfa50(1);
PGA90vec84 = PGA90vec*alfa84(1);
PGV90vec50 = PGV90vec*alfa50(2);
PGV90vec84 = PGV90vec*alfa84(2);
PGD90vec50 = PGD90vec*alfa50(3);
PGD90vec84 = PGD90vec*alfa84(3);

%COMPONENTE 360º
PGA360vec50 = PGA360vec*alfa50(1);
PGA360vec84 = PGA360vec*alfa84(1);
PGV360vec50 = PGV360vec*alfa50(2);
PGV360vec84 = PGV360vec*alfa84(2);
PGD360vec50 = PGD360vec*alfa50(3);
PGD360vec84 = PGD360vec*alfa84(3);

%COMPONENTE UP
PGAUPvec50 = PGAUPvec*alfa50(1);
PGAUPvec84 = PGAUPvec*alfa84(1);
PGVUPvec50 = PGVUPvec*alfa50(2);
PGVUPvec84 = PGVUPvec*alfa84(2);
PGDUPvec50 = PGDUPvec*alfa50(3);
PGDUPvec84 = PGDUPvec*alfa84(3);

%%%%%%%%%%%%%%%%% LINEA BASE MEDIANA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LB90_50 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 
LB360_50 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 
LBUP_50 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 

for i =1:length(Ttetra)
   LB90_50(i) = min([PGA90vec50(i) PGV90vec50(i) PGD90vec50(i)]);
end

for i =1:length(Ttetra)
   LB360_50(i) = min([PGA360vec50(i) PGV360vec50(i) PGD360vec50(i)]);
end

for i =1:length(Ttetra)
   LBUP_50(i) = min([PGAUPvec50(i) PGVUPvec50(i) PGDUPvec50(i)]);
end

%%%%%%%%%%%%%%%%% LINEA BASE PERCENTIL 84 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LB90_84 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 
LB360_84 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 
LBUP_84 = zeros(length(Ttetra),1); %% vector para la construccion del la linea base componente 90 

for i =1:length(Ttetra)
   LB90_84(i) = min([PGA90vec84(i) PGV90vec84(i) PGD90vec84(i)]);
end

for i =1:length(Ttetra)
   LB360_84(i) = min([PGA360vec84(i) PGV360vec84(i) PGD360vec84(i)]);
end

for i =1:length(Ttetra)
   LBUP_84(i) = min([PGAUPvec84(i) PGVUPvec84(i) PGDUPvec84(i)]);
end



%% ploteo Las Pseudo Velocidades en formato tetralogaritmico


 
figure(2)
subplot(1,3,1)
Representacion_tetralogaritmica()
title('Componente 90º')
hold on 
loglog(Ttetra,Sv(:,1),'b','lineWidth',2);
loglog(Ttetra,LB90,'c','lineWidth',2);
loglog(Ttetra,LB90_50,'r','lineWidth',2);
loglog(Ttetra,LB90_84,'g','lineWidth',2);
legend('Espectro Pseudo-velocidades (Azul)','Linea base (Cyan)','Linea base - mediana (Rojo)','Linea base - percentil 84% (Verde)');



subplot(1,3,2)
Representacion_tetralogaritmica()
title('Componente 360º')
hold on 
loglog(Ttetra,Sv(:,2),'r','linewidth',2);
loglog(Ttetra,LB360,'m','lineWidth',2);
loglog(Ttetra,LB360_50,'g','lineWidth',2);
loglog(Ttetra,LB360_84,'b','lineWidth',2);
legend('Espectro Pseudo-velocidades (Rojo)','Linea base (Magenta)','Linea base - mediana (Verde)','Linea base - percentil 84% (Azul)');



subplot(1,3,3)
Representacion_tetralogaritmica()
title('Componente up')
hold on 
loglog(Ttetra,Sv(:,3),'g','linewidth',2);
loglog(Ttetra,LBUP,'y','lineWidth',2);
loglog(Ttetra,LBUP_50,'r','lineWidth',2);
loglog(Ttetra,LBUP_84,'b','lineWidth',2);
legend('Espectro Pseudo-velocidades (Verde)','Linea base (Amarillo)','Linea base - mediana (Rojo)','Linea base - percentil 84% (Azul)');




