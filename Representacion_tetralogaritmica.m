function Representacion_tetralogaritmica()
%************************************************************************
%Rutina para generar representación tri-partita de espectros de respuesta 
%Es necesario leer de archivo de espectro de respuesta
%************************************************************************
for k=.00001:.00001:.0001
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=.0001:.0001:.001
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=.001:.001:.01
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
xlabel(' Periodo (s)')
ylabel(' Velocidad Espectral (cm/s)')
for k=.01:.01:.1
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=.1:.1:1
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=1:1:10
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=10:10:100
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=100:100:1000
x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
hold on
end
for k=1000:1000:10000
    x=0.01:1:100;
t=log(2*pi*k)-log(x);
y=exp(t);
loglog(x,y,'k'),grid on
hold on
t=log(k*9.81/(2*pi))+log(x);
y=exp(t);
loglog(x,y,'k')
end
axis([0.01 10 0.02 500])
% d=xlsread('svdata');
% sv='sv.out'
% d=load(sv)
% plot(d(:,1),100*d(:,2),'k')
% plot(d(:,1),100*d(:,3),'k')
% plot(d(:,1),100*d(:,4),'k')
% plot(d(:,1),100*d(:,5),'k')
% plot(d(:,1),100*d(:,6),'k')
% text(0.2,0.02,'0.001');
% text(0.6,0.1,'0.01');
% text(2,0.3,'0.1');
% text(7,1,'1');
% text(20,3,'10');
% text(80,10,'100')
% text(20,1,'Sd  (cm)')
% xlabel(' periodo (s)')
% ylabel(' Sv  (cm/s)')
% text(0.01,200,'100')
% text(0.01,20,'10')
% text(0.01,2,'1')
% text(0.02,0.4,'0.1')
% text(0.07,0.1,'0.01')
% text(.02,0.8,'sa/g')
% gtext(' sin amortiguamiento')
% gtext(' amortiguamiento=2%')
% gtext(' amortiguamiento=5%')
% gtext(' amortiguamiento=10%')
% gtext(' amortiguamiento=20%')
end
