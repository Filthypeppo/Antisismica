function [PA,PV,PD,tt,a,v,x,Q]=Beta_Newmark_NoLineal(zi,B,Wo,xo,vo,t,Sg,ko,Qy)
% PROYECTO DE INVESTIGACION POTENCIAL DESTRUCTIVO GENERADO POR EL TERREMOTO
%              DEL 16 DE ABRIL DE 2016 - PRESENTADO AL IPGH
%                      ARGENTINA - ECUADOR - CHILE
%                            Abril 2019        

% Por Abel Zambrano, ULEAM
% Revisado por Juan Vielma, PUCV
% function [PA,PV,PD,tt,a,v,x,Q]=Beta_Newmark_NoLineal(zi,B,Wo,xo,vo,t,Sg,ko,Qy)
% La función calcula la respuesta no lineal paso a paso por el método Beta 
% de Newmark para sistemas de un grado de libertad considerando un sistema 
% elastoplástico

% Datos de entrada:
% zi=razón de amortiguamiento
% B=Constante Beta entre 1/8-1/4
% Wo=Frecuencia natural del sistema
% xo=Desplazamiento inicial
% vo=Velocidad inicial
% t=tiempo del acelerograma
% Sg=Aceleración del acelerograma
% ko=Rigidez inicial del sistema [T/cm]
% Qy=Rigidez de fluencia [T]

% Datos de salida:
% PA=Aceleración máxima del sistema elastoplástico
% PVe=Velocidad máxima del sistema elastoplástico
% PDe=Desplazamiento máximo del sistema elastoplástico
% t=Intantes de tiempo para las respuestas de acel.,vel.,desp.
% a=Respuesta paso a paso de aceleración del sistema elastoplástico
% v=Respuesta paso a paso de velocidad sistema elastoplástico
% x=Respuesta paso a paso de desplazamiento sistema elástico asociado
% Q=Respuesta paso a paso de la resistencia del sistema elástico asociado


%% CÁLCULO DE LA RESPUESTA PASO A PASO EN LA ZONA ELÁSTICA INICIAL

dt=t(2)-t(1);
n=length(t);

% Primer valor de aceleración
ao=-Sg(1)-2*zi*Wo*vo-Wo^2*xo;

a=ao;
v=vo;
x=xo;
Q=ko*xo;

% Desplazamiento de fluencia
xy=Qy/ko;
tt=t(1);
Sgt=Sg(1);

aux1=0;
for k=2:n
i=k+aux1;
a(i)=((1+zi*Wo*dt+B*Wo^2*dt^2)^-1)*(-Sg(k)-Wo^2*x(i-1)-(Wo^2*dt+2*zi*Wo)*v(i-1)-((0.5-B)*Wo^2*dt^2+zi*Wo*dt)*a(i-1));
v(i)=v(i-1)+dt/2*(a(i-1)+a(i));
x(i)=x(i-1)+v(i-1)*dt+((0.5-B)*a(i-1)+B*a(i))*dt^2;

tt(i)=t(k);
Sgt(i)=Sg(k);

if abs(x(i))<xy
Q(i)=ko*x(i);
else
    
% CÁLCULO DEL INSTANTE DONDE Q=Qy

t1=t(k-1); t2=t(k);
x1=x(i-1); x2=x(i); % CORRECCION x2=v(i)
v1=v(i-1); a1=a(i-1);
Sg1=Sg(k-1); Sg2=Sg(k);

if x2>0;
   xm=xy;
else
   xm=-xy;
end

to=t1-(t2-t1)*(xm-x1)/(x2-x1); % Interpolación simple para calcular el tiempo donde Q=Qy (v. preliminar)
dto=to-t1;
Sgo=(Sg2-Sg1)/(t2-t1)*dto+Sg1; % Valor de aceleración Sg en el tiempo donde Q=Qy (v. preliminar)

% Aproximación del valor de aceleración para el instante donde Q=Qy
ao=-Sgo-2*zi*Wo*v1-Wo^2*x1; % Se utiliza la fórmula iterativa
vo=v1+dto/2*(a1+ao);
xo=x1+v1*dto+((0.5-B)*a1+B*ao)*dto^2;
ea=0.0001; % Tolerancia del método iterativo
maxite=100; % Máxima iteración

% Cálculo del instante donde Q=Qy aplicando el método iterativo
for j=2:maxite
to(j)=t1+(to(j-1)-t1)*(xm-x1)/(xo(j-1)-x1);
dto(j)=to(j)-t1;
Sgo(j)=(Sg2-Sg1)/(t2-t1)*dto(j)+Sg1; % Se interpola el valor de aceleración del acelerograma
ao(j)=-Sgo(j)-2*zi*Wo*vo(j-1)-Wo^2*xo(j-1);
vo(j)=v1+dto(j)/2*(a1+ao(j));
xo(j)=x1+v1*dto(j)+((0.5-B)*a1+B*ao(j))*dto(j)^2;

if abs(ao(j)-ao(j-1))<=ea
    break
end

end

%fprintf ('Instante donde Q=Qy')
t1+dto(end);

tt(i)=t1+dto(end);
Sgt(i)=Sgo(end);

% Cálculo de la respuesta en el tiempo donde Q=Qy
dt1=t2-t1-dto(end);

a(i)=ao(end);
v(i)=vo(end);
x(i)=xo(end);
Q(i)=ko*x(i);

break
end

end


%% CÁLCULO DE LA RESPUESTA PASO A PASO EN ZONA PLÁSTICA Y ELÁSTICA DESFASADA

% Inicio de la zona plástica (continuación del análisis)
ini=i;

%zona=1; % Zona plástica
%zona=2; % Zona elástica desfasada

% Cálculo automático desde el tercer paso en adelante
if v(ini-1)>0
    aux(ini-1)=1;
    else
    aux(ini-1)=2;
end

xu=[]; % Aquí se guardarán los desplazamientos últimos
ub_xu=[];
aux2=0; % Auxiliar para añadir los tiempos de velocidad nula y cuando se alcanza la fluencia

zona=1; % Se inicia en zona 1, puesto que se continúa el análisis con la zona plástica y luego elástica desfasada

%% ANÁLISIS EN LA ZONA PLÁSTICA Y ELÁSTICA DESFASADA

for k=ini:n

%==========================================================================
% ZONA PLÁSTICA
%==========================================================================

if zona==1
    
i=1+k+aux2; % Se suma 1 para añadir el tiempo adicional de la primera fluencia

tt(i)=t(k); % Guarda el vector de tiempo considerando los tiempos adicionales de Q=Qy y v=0
Sgt(i)=Sg(k);

dt=tt(i)-tt(i-1); % Calcula del diferencial de tiempo
    
Qm=Q(i-1);

a(i)=((1+zi*Wo*dt)^-1)*(-2*zi*Wo*v(i-1)-zi*Wo*dt*a(i-1)-Wo^2*Q(i-1)/Qy-Sg(k));
v(i)=v(i-1)+dt/2*(a(i-1)+a(i));


% Guarda un código si la velocidad es positiva aux=1, caso contratio aux=2
if v(i)>0
    aux(k)=1;
else
    aux(k)=2;
end

% Verifica si la velocidad no cambia de signo
if aux(k)==aux(k-1)

x(i)=x(i-1)+v(i-1)*dt+((0.5-B)*a(i-1)+B*a(i))*dt^2;

Q(i)=ko*Qm/ko;

else

% Cálculo del instante donde la velocidad es nula
t1=t(k-1); t2=t(k);
v1=v(i-1); v2=v(i);
a1=a(i-1); x1=x(i-1);
Sg1=Sg(k-1); Sg2=Sg(k);

to=t1-v1*(t2-t1)/(v2-v1); % Interpolación simple para calcular el tiempo de velocidad nula (v. preliminar)
dto=to-t1;
Sgo=(Sg2-Sg1)/(t2-t1)*dto+Sg1; % Valor de aceleración Sg en el tiempo donde la velocidad es nula (v. preliminar)

% Aproximación del valor de aceleración para el instante de velocidad nula
if Qm>0
   xm=xy;
else
   xm=-xy;
end
ao=((1+zi*Wo*dto)^-1)*(-2*zi*Wo*v1-zi*Wo*dto*a1-Wo^2*xm-Sgo);

vo=v1+dto/2*(a1+ao);
xo=x1+v1*dto+((0.5-B)*a1+B*ao)*dto^2;
ea=0.0001; % Tolerancia del método iterativo
maxite=100; % Máxima iteración

% Cálculo del instante donde la velocidad es nula aplicando el método
% iterativo
for j=2:maxite
dto(j)=-2*v1/(a1+ao(j-1)); % Se despeja dt de la fórmula de la velocidad haciendo v(n+1)=0
Sgo(j)=(Sg2-Sg1)/(t2-t1)*dto(j)+Sg1; % Se interpola el valor de aceleración del acelerograma
ao(j)=((1+zi*Wo*dto(j))^-1)*(-2*zi*Wo*v1-zi*Wo*dto(j)*a1-Wo^2*xm-Sgo(j)); % Se utiliza la fórmula la rango plástico
vo(j)=v1+dto(j)/2*(a1+ao(j)); % Se utiliza la fórmula de velocidad y desplazamiento convencional
xo(j)=x1+v1*dto(j)+((0.5-B)*a1+B*ao(j))*dto(j)^2;

if abs(ao(j)-ao(j-1))<=ea
    break
end
end

tt(i)=t1+dto(end);
Sgt(i)=Sgo(end);

%Cálculo de la respuesta en el tiempo t2
dt1=t2-t1-dto(end);

a(i)=ao(end);
v(i)=vo(end);
x(i)=xo(end);
Q(i)=ko*xm;

% Guarda el valor del desplazamiento último
xu(end+1)=x(i);
ub_xu(end+1)=i;

tt(i+1)=t2;
Sgt(i+1)=Sg2;

if Qm<0;
    xr=xo(end)-(xu(end)+xy);
else
    xr=xo(end)-(xu(end)-xy);
end

a(i+1)=((1+zi*Wo*dt1+B*Wo^2*dt1^2)^-1)*(-Sg2-Wo^2*xr-(Wo^2*dt1+2*zi*Wo)*vo(end)-((0.5-B)*Wo^2*dt1^2+zi*Wo*dt1)*ao(end));
v(i+1)=vo(end)+dt1/2*(ao(end)+a(i+1));
x(i+1)=xo(end)+vo(end)*dt1+((0.5-B)*ao(end)+B*a(i+1))*dt1^2;

if Qm<0;
    xr=x(end)-(xu(end)+xy);
else
    xr=x(end)-(xu(end)-xy);
end

Q(i+1)=ko*xr;

aux2=aux2+1;

if v(i+1)>0
    aux(k-1)=1;
    else
    aux(k-1)=2;
end


% Pasa de la zona plástica a la zona elástica desfasada
zona=2;

end

%==========================================================================
% ZONA ELÁSTICA DESFASADA
%========================================================================== 
elseif zona==2 
    
i=1+k+aux2; % Se suma por el tiempo adicional de la primera fluencia

tt(i)=t(k);
Sgt(i)=Sg(k);

dt=tt(i)-tt(i-1);

if Qm<0
    xr=x(i-1)-(xu(end)+xy);
else
    xr=x(i-1)-(xu(end)-xy);
end

tt(i)=t(k);
Sgt(i)=Sg(k);

a(i)=((1+zi*Wo*dt+B*Wo^2*dt^2)^-1)*(-Sg(k)-Wo^2*xr-(Wo^2*dt+2*zi*Wo)*v(i-1)-((0.5-B)*Wo^2*dt^2+zi*Wo*dt)*a(i-1));
v(i)=v(i-1)+dt/2*(a(i-1)+a(i));
x(i)=x(i-1)+v(i-1)*dt+((0.5-B)*a(i-1)+B*a(i))*dt^2;

if Qm<0
    xr=x(i)-(xu(end)+xy);
else
    xr=x(i)-(xu(end)-xy);
end

Q(i)=ko*xr;

if abs(Q(i))>Qy

% CÁLCULO DEL INSTANTE DONDE Q=Qy
t1=t(k-1); t2=t(k);
Q1=Q(i-1); Q2=Q(i);
v1=v(i-1); a1=a(i-1); x1=x(i-1);
Sg1=Sg(k-1); Sg2=Sg(k);

if Q2>0
   Qm1=Qy; 
else
   Qm1=-Qy;
end

to=t1+(t2-t1)*(Qm1-Q1)/(Q2-Q1); % Interpolación simple para calcular el tiempo donde Q=Qy (v. preliminar)
dto=to-t1;
Sgo=(Sg2-Sg1)/(t2-t1)*dto+Sg1; % Valor de aceleración Sg en el tiempo donde Q=Qy (v. preliminar)

% Aproximación del valor de aceleración para el instante donde Q=Qy

ao=((1+zi*Wo*dto+B*Wo^2*dto^2)^-1)*(-Sgo-Wo^2*(x1-(xu(end)-Qm/ko))-(Wo^2*dto+2*zi*Wo)*v1-((0.5-B)*Wo^2*dto^2+zi*Wo*dto)*a1);

%ao=-Sgo-2*zi*Wo*v1-Wo^2*(x1-(xu(end)-Qm/ko)) % Se utiliza la fórmula iterativa

vo=v1+dto/2*(a1+ao);
xo=x1+v1*dto+((0.5-B)*a1+B*ao)*dto^2;
Qo=ko*(xo-(xu(end)-Qm/ko));

ea=0.00001; % Tolerancia del método iterativo
maxite=100; % Máxima iteración

% Cálculo del instante donde Q=Qy aplicando el método iterativo
for j=2:maxite
if abs(Qo(j-1))>abs(Qm1)
    t1p=t1; t2p=to(end);
    Q1p=Q1; Q2p=Qo(end);
to(j)=t1p+(t2p-t1p)*(Qm1-Q1p)/(Q2p-Q1p);
dto(j)=to(j)-t1;
Sgo(j)=(Sg2-Sg1)/(t2-t1)*dto(j)+Sg1; % Se interpola el valor de aceleración del acelerograma
%ao(j)=-Sgo(j)-2*zi*Wo*vo(j-1)-Wo^2*(xo(j-1)-(xu(end)-Qm/ko)); %Iterativo
ao(j)=((1+zi*Wo*dto(j)+B*Wo^2*dto(j)^2)^-1)*(-Sgo(j)-Wo^2*(x1-(xu(end)-Qm/ko))-(Wo^2*dto(j)+2*zi*Wo)*v1-((0.5-B)*Wo^2*dto(j)^2+zi*Wo*dto(j))*a1);
vo(j)=v1+dto(j)/2*(a1+ao(j));
xo(j)=x1+v1*dto(j)+((0.5-B)*a1+B*ao(j))*dto(j)^2;
Qo(j)=ko*(xo(j)-(xu(end)-Qm/ko)); % Cambio signo Qm
if abs(Qo(j)-Qo(j-1))<=ea
    break
end
else
t1p=to(end); t2p=t2;
Q1p=Qo(end); Q2p=Q2;
to(j)=t1p+(t2p-t1p)*(Qm1-Q1p)/(Q2p-Q1p);
dto(j)=to(j)-t1;
Sgo(j)=(Sg2-Sg1)/(t2-t1)*dto(j)+Sg1; % Se interpola el valor de aceleración del acelerograma
%ao(j)=-Sgo(j)-2*zi*Wo*vo(j-1)-Wo^2*(xo(j-1)-(xu(end)-Qm/ko)); %Iterativo
ao(j)=((1+zi*Wo*dto(j)+B*Wo^2*dto(j)^2)^-1)*(-Sgo(j)-Wo^2*(x1-(xu(end)-Qm/ko))-(Wo^2*dto(j)+2*zi*Wo)*v1-((0.5-B)*Wo^2*dto(j)^2+zi*Wo*dto(j))*a1);
vo(j)=v1+dto(j)/2*(a1+ao(j));
xo(j)=x1+v1*dto(j)+((0.5-B)*a1+B*ao(j))*dto(j)^2;
Qo(j)=ko*(xo(j)-(xu(end)-Qm/ko)); % Cambio signo Qm 
if abs(Qo(j)-Qo(j-1))<=ea
    break
end
end
end

tt(i)=t1+dto(end);
Sgt(i)=Sgo(end);

% Cálculo de la respuesta en el tiempo donde Q=Qy
dt1=t2-t1-dto(end);

a(i)=ao(end);
v(i)=vo(end);
x(i)=xo(end);
Q(i)=ko*(x(i)-(xu(end)-Qm/ko)); % Cambio de signo

% Cálculo en el tiempo t2
tt(i+1)=t2;
Sgt(i+1)=Sg2;
a(i+1)=((1+zi*Wo*dt1)^-1)*(-2*zi*Wo*v(i)-zi*Wo*dt1*a(i)-Wo^2*Q(i)/ko-Sg(k));
v(i+1)=v(i)+dt1/2*(a(i)+a(i+1));
x(i+1)=x(i)+v(i)*dt1+((0.5-B)*a(i)+B*a(i+1))*dt1^2;
Q(i+1)=ko*Q(i)/ko;

Qm=Q(i+1);

aux2=aux2+1; % Cuenta la primera fluencia    

%Grabacion del signo de la velocidad

if v(i+1)>0
    aux(k)=1;
    else
    aux(k)=2;
end

% De la zona elástica desfasada pasa a la zona plástica
zona=1;

%==========================================================================
% CASO PARTICULAR.- Donde se pasa del rango plástico al elástico desfasado
% en un corto intervalo de tiempo.
if v(i)>0
    aux(k-1)=1;
else
    aux(k-1)=2;
end

if aux(k-1)~=aux(k)

% Cálculo del instante donde la velocidad es nula en un rango plástico
% corto

t1=tt(i); t2=tt(i+1);
v1=v(i); v2=v(i+1);
a1=a(i); x1=x(i);
Sg1=Sgt(i); Sg2=Sgt(i+1);

to=t1-v1*(t2-t1)/(v2-v1); % Interpolación simple para calcular el tiempo de velocidad nula (v. preliminar)
dto=to-t1;
Sgo=(Sg2-Sg1)/(t2-t1)*dto+Sg1; % Valor de aceleración Sg en el tiempo donde la velocidad es nula (v. preliminar)

% Aproximación del valor de aceleración para el instante de velocidad nula
if Q(i)>0
    Qm=Qy;
else
    Qm=-Qy;
end

ao=((1+zi*Wo*dto)^-1)*(-2*zi*Wo*v1-zi*Wo*dto*a1-Wo^2*Qm/ko-Sgo);
vo=v1+dto/2*(a1+ao);
xo=x1+v1*dto+((0.5-B)*a1+B*ao)*dto^2;
ea=0.0001; % Tolerancia del método iterativo
maxite=100; % Máxima iteración

% Cálculo del instante donde la velocidad es nula aplicando el método
% iterativo
for j=2:maxite
dto(j)=-2*v1/(a1+ao(j-1)); % Se despeja dt de la fórmula de la velocidad haciendo v(n+1)=0
Sgo(j)=(Sg2-Sg1)/(t2-t1)*dto(j)+Sg1; % Se interpola el valor de aceleración del acelerograma
ao(j)=((1+zi*Wo*dto(j))^-1)*(-2*zi*Wo*v1-zi*Wo*dto(j)*a1-Wo^2*Qm/ko-Sgo(j)); % Se utiliza la fórmula la rango plástico
vo(j)=v1+dto(j)/2*(a1+ao(j)); % Se utiliza la fórmula de velocidad y desplazamiento convencional
xo(j)=x1+v1*dto(j)+((0.5-B)*a1+B*ao(j))*dto(j)^2;

if abs(ao(j)-ao(j-1))<=ea
    break
end
end


tt(i+1)=t1+dto(end);
Sgt(i+1)=Sgo(end);

%Cálculo de la respuesta en el tiempo t2
dt1=t2-t1-dto(end);

a(i+1)=ao(end);
v(i+1)=vo(end);
x(i+1)=xo(end);
Q(i+1)=ko*Qm/ko;

% Guarda el valor del desplazamiento último
xu(end+1)=x(i+1);
ub_xu(end+1)=i+1;

% Las respuestas en el nuevo tiempo "t2" t(i+2)
tt(i+2)=t2;
Sgt(i+2)=Sg2;

if Qm<0;
    xr=xo(end)-(xu(end)+xy);
else
    xr=xo(end)-(xu(end)-xy);
end

a(i+2)=((1+zi*Wo*dt1+B*Wo^2*dt1^2)^-1)*(-Sg2-Wo^2*xr-(Wo^2*dt1+2*zi*Wo)*vo(end)-((0.5-B)*Wo^2*dt1^2+zi*Wo*dt1)*ao(end));
v(i+2)=vo(end)+dt1/2*(ao(end)+a(i+1));
x(i+2)=xo(end)+vo(end)*dt1+((0.5-B)*ao(end)+B*a(i+1))*dt1^2;

if Qm<0;
    xr=x(end)-(xu(end)+xy);
else
    xr=x(end)-(xu(end)-xy);
end

Q(i+2)=ko*xr;

aux2=aux2+1;

% Esta parte me surge duda
if v(i+2)>0
    aux(k-1)=1;
    else
    aux(k-1)=2;
end

zona=2;

end

end

%==========================================================================

end

end

figure1=figure;
plot(x,Q); 
xlabel({'Displacement [cm]'});
ylabel({'Strength [Ton]'});
title('Elasto-plastic Behavior');

%% RESPUESTAS MÁXIMAS
PA=max(max(a),abs(min(a)));
PV=max(max(v),abs(min(v)));
PD=max(max(x),abs(min(x)));

%fprintf('       RESUMEN DE RESULTADOS \n');
fprintf('\n Aceleración máxima no lineal [cm/s2] \n %3.4f\n',PA);
fprintf('\n Velocidad máxima no lineal [cm/s] \n %3.4f\n',PV);
fprintf('\n Desplazamiento máximo no lineal [cm] \n %3.4f\n',PD);

%[tt',Sgt',a',v',x',Q']

end


