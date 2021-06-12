%% PROYECTO DE INVESTIGACION POTENCIAL DESTRUCTIVO GENERADO POR EL TERREMOTO
%              DEL 16 DE ABRIL DE 2016 - PRESENTADO AL IPGH
%                      ARGENTINA - ECUADOR - CHILE

%           SISTEMAS NO LINEALES DE UN GRADO DE LIBERTAD
%       RESPUESTA PASO A PASO USANDO EL M�TODO BETA DE NEWMARK
%                           Abril 2019

% Por Abel Zambrano, ULEAM
% Revisado por Juan Vielma, PUCV
% Inicializaci�n
clc
clear all

%% DATOS DE ENTRADA:
zi=0.05; % Raz�n de amortiguamiento 
B=0.25; % Constante Beta (Entre 1/6 a 1/4)
Wo=12.56637061; % Frecuencia natural [rad/s]
ko=30; % Rigidez inicial del sistema [T/cm]
Qy=30; % Rigidez de fluencia [T]
% Condiciones iniciales
xo=0; % Desplazamiento inicial
vo=0; % Velocidad inicial

%% CARGAR EL ACELEROGRAMA
% Unidades del registro [t(s),Sg(cm/s2)]
load 'Comp360.txt'
Acelerograma=Comp360;
t=Acelerograma(:,1);
Sg=Acelerograma(:,2)*981;

%% AN�LISIS DEL SISTEMA EL�STICO ASOCIDADO
% Llama a la rutina "Beta_Newmark" para calcular la respuesta lineal en el
% tiempo con el m�todo de Beta de Newmark
[PAe,PVe,PDe,te,ae,ve,xe]=Beta_Newmark(zi,B,Wo,xo,vo,t,Sg);

%% AN�LISIS DEL SISTEMA ELASTOPL�STICO
% Llama a la rutina "Beta_Newmark_NoLineal" para calcular respuesta no lineal en el
% tiempo del sistema el�stopl�stico con el m�todo de Beta de Newmark
[PA,PV,PD,tt,a,v,x,Q]=Beta_Newmark_NoLineal(zi,B,Wo,xo,vo,t,Sg,ko,Qy);

%% RESPUESTAS
miu=PD*ko/Qy;
fprintf('\n Factor de ductilidad \n %3.2f\n',miu);
R=ko*PDe/Qy;
fprintf('\n Factor de reducci�n de resistencia \n %3.2f\n',R);


%% Create figure displacement response
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'on');

% Create plot
plot(te,xe,'DisplayName','Elastic','Color',[0 0 1]); %b


% Create plot
plot(tt,x,'DisplayName','Nonlinear','Color',[1 0 0]); %r

% Create xlabel
xlabel({'Time [sec]'});

% Create ylabel
ylabel({'Displacement [cm]'});

% Create title
title('Elastic and nonlinear response');

% Create legend
legend(axes1,'show');

