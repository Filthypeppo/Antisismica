function [PAe,PVe,PDe,te,ae,ve,xe]=Beta_Newmark(zi,B,Wo,xo,vo,t,Sg)
% PROYECTO DE INVESTIGACION POTENCIAL DESTRUCTIVO GENERADO POR EL TERREMOTO
%              DEL 16 DE ABRIL DE 2016 - PRESENTADO AL IPGH
%                      ARGENTINA - ECUADOR - CHILE
%                           Marzo 2019

% Por Abel Zambrano, ULEAM
% Revisado por Juan Vielma, PUCV
% function [PA,PV,PD,ae,ve,xe]=Beta_Newmark(zi,B,Wo,xo,vo,t,Sg)
% La función calcula la respuesta lineal paso a paso por el método Beta de 
% Newmark para sistemas de un grado de libertad 

% Datos de entrada:
% zi=razón de amortiguamiento
% B=Constante Beta entre 1/8-1/4
% Wo=Frecuencia natural del sistema
% xo=Desplazamiento inicial
% vo=Velocidad inicial
% t=tiempo del acelerograma
% Sg=Aceleración del acelerograma

% Datos de salida:
% PAe=Aceleración máxima elástica
% PVe=Velocidad máxima elástica
% PDe=Desplazamiento máximo elástico
% te=Intantes de tiempo para las respuestas de acel.,vel.,desp.
% ae=Respuesta pasao a paso de aceleración sistema elástico asociado
% ve=Respuesta de velocidad sistema elástico asociado
% xe=Respuesta de desplazamiento sistema elástico asociado

%--------------------------------------------------------------------------

dt=t(2)-t(1);
n=length(t);

% CÁLCULO DE LA RESPUESTA PASO A PASO
% Primer valor de aceleración
ao=-Sg(1)-2*zi*Wo*vo-Wo^2*xo;

a=ao;
v=vo;
x=xo;

% Segundo cálculo manual
a(2)=((1+zi*Wo*dt+B*Wo^2*dt^2)^-1)*(-Sg(2)-Wo^2*x(1)-(Wo^2*dt+2*zi*Wo)*v(1)-((0.5-B)*Wo^2*dt^2+zi*Wo*dt)*a(1));
v(2)=v(1)+dt/2*(a(1)+a(2));
x(2)=x(1)+v(1)*dt+((0.5-B)*a(1)+B*a(2))*dt^2;
tt=t(1:2);

% Cálculo automático desde el tercer paso en adelante
if v(2)>0
    aux(2)=1;
    else
    aux(2)=2;
end

aux2=0;
for k=3:n
i=k+aux2;
a(i)=((1+zi*Wo*dt+B*Wo^2*dt^2)^-1)*(-Sg(k)-Wo^2*x(i-1)-(Wo^2*dt+2*zi*Wo)*v(i-1)-((0.5-B)*Wo^2*dt^2+zi*Wo*dt)*a(i-1));
v(i)=v(i-1)+dt/2*(a(i-1)+a(i));
tt(i)=t(k);
if v(i)>0
    aux(k)=1;
    else
    aux(k)=2;
end

if aux(k)==aux(k-1)
x(i)=x(i-1)+v(i-1)*dt+((0.5-B)*a(i-1)+B*a(i))*dt^2;
else
    
% CÁLCULO DEL INSTANTE DONDE LA VELOCIDAD ES NULA

t1=t(k-1); t2=t(k);
v1=v(i-1); v2=v(i);
a1=a(i-1); x1=x(i-1);
Sg1=Sg(k-1); Sg2=Sg(k);

to=t1-v1*(t2-t1)/(v2-v1); % Interpolación simple para calcular el tiempo de velocidad nula (v. preliminar)
dto=to-t1;
Sgo=(Sg2-Sg1)/(t2-t1)*dto+Sg1; % Valor de aceleración Sg en el tiempo donde la velocidad es nula (v. preliminar)

% Aproximación del valor de aceleración para el instante de velocidad nula
ao=-Sgo-Wo^2*x1; % Se utiliza la fórmula iterativa con velocidad nula
vo=v1+dto/2*(a1+ao);
xo=x1+v1*dto+((0.5-B)*a1+B*ao)*dto^2;
ea=0.0001; % Tolerancia del método iterativo
maxite=20; % Máxima iteración

% Cálculo del instante donde la velocidad es nula aplicando el método
% iterativo
for j=2:maxite
dto(j)=-2*v1/(a1+ao(j-1)); % Se despeja dt de la fórmula de la velocidad haciendo v(n+1)=0
Sgo(j)=(Sg2-Sg1)/(t2-t1)*dto(j)+Sg1; % Se interpola el valor de aceleración del acelerograma
ao(j)=-Sgo(j)-Wo^2*xo(j-1); % Se utiliza la fórmula de la aceleración iterativa sabiendo que la v(n+1)=0
vo(j)=v1+dto(j)/2*(a1+ao(j)); % Se utiliza la fórmula de velocidad y desplazamiento convencional
xo(j)=x1+v1*dto(j)+((0.5-B)*a1+B*ao(j))*dto(j)^2;

if abs(ao(j)-ao(j-1))<=ea
    break
end

end


%fprintf ('Instante donde la velocidad es nula')
%t1+dto(end)

tt(i)=t1+dto(end);

% Cálculo de la respuesta en el tiempo t2
dt1=t2-t1-dto(end);

a(i)=ao(end);
v(i)=vo(end);
x(i)=xo(end);

tt(i+1)=t2;
a(i+1)=((1+zi*Wo*dt1+B*Wo^2*dt1^2)^-1)*(-Sg2-Wo^2*xo(end)-(Wo^2*dt1+2*zi*Wo)*vo(end)-((0.5-B)*Wo^2*dt1^2+zi*Wo*dt1)*ao(end));
v(i+1)=vo(end)+dt1/2*(ao(end)+a(i+1));
x(i+1)=xo(end)+vo(end)*dt1+((0.5-B)*ao(end)+B*a(i+1))*dt1^2;

aux2=aux2+1;

end

end

% Guarda los resultados
te=tt; % Tiempo para los intantes de tiempo para las respuestas de a,v,x
ae=a; % Respuesta de aceleración sistema elástico asociado
ve=v; % Respuesta de velocidad sistema elástico asociado
xe=x; % Respuesta de desplazamiento sistema elástico asociado


% Create figure displacement response
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'on');

% Create plot
subplot(1,3,1),plot(tt,a,'DisplayName','Acceleration','Color',[1 1 0]); %y[1 1 0]  b[0 0 1]
title('Aceleracion')
xlabel({'Time [sec]'})
ylabel({'Elastic Response [cm,sec]'});

% Create plot
subplot(1,3,2),plot(tt,v,'DisplayName','Velocity','Color',[0 1 0]); %g
title('Velocidad')
xlabel({'Time [sec]'})
ylabel({'Elastic Response [cm,sec]'});

% Create plot
subplot(1,3,3),plot(tt,x,'DisplayName','Displacement','Color',[1 0 0]); %r
title('Desplazamiento')
xlabel({'Time [sec]'})
ylabel({'Elastic Response [cm,sec]'});

% Create xlabel
%xlabel({'Time [sec]'});

% Create ylabel
%ylabel({'Elastic Response [cm,sec]'});

% Create title
%title('Elastic Response');

% Create legend
%legend(axes1,'show');


PAe=max(max(a),abs(min(a)));
PVe=max(max(v),abs(min(v)));
PDe=max(max(x),abs(min(x)));

fprintf('       RESUMEN DE RESULTADOS \n');
fprintf('\n Aceleración máxima elástica [cm/s2] \n %3.4f\n',PAe);
fprintf('\n Velocidad máxima elástica [cm/s] \n %3.4f\n',PVe);
fprintf('\n Desplazamiento máximo elástico [cm] \n %3.4f\n',PDe);

end
