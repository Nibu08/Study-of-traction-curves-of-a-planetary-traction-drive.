%% Method for estimating traction curves - creep theory.

%% 
clc
clear all
Cr=0.04;
V_oil=linspace(2,9.2,100);
R1=140e-3%input('Rotational radius of drive side rolling element:');
R2=120e-3%input('Rotational radius of driven side rolling element:');
a=R1/8%input('Long radius of contact ellipse :');
b=R1/2%input('short radius of contact ellipse :');
U1=10%input('Peripheral Velocity of drive side rolling element : ');
w1=6%input('angular Velocity of drive side rolling element : ');
w2=3%input('angular Velocity of driven side rolling element : ');
U2=15%input('Peripheral Velocity of driven side rolling element : ');
Fc=20e3 %input('Clamping Force (Axial Force) : ' );
th1=45%input('Angle between tangent plane &Rotate axis of drive side rolling element :');
th2=36%input('Angle between tangent plane &Rotate axis of driven side rolling element :');
P_max =3e9%input('Max contact pressure  :');%   Max contact pressure
G=150e6%input(' Effictive Shear Modules :');
v=0.3%input('poison ratio of rolling element');
x = linspace(0,a/2,100); 
y = linspace(0,b/2,100); 
[X,Y] = meshgrid(x,y);
tau_e=[]
%% Peripheral Velocity of drive side rolling element
wsp= (w1.*sind(th1))-(w2.*sind(th2));% Spin Angular velocity 
ux=(U1.*cosd(th1))-(U2.*cosd(th2));% Resulant velocity  at x axis 
uy=(U1.*sind(th1))-(U2.*sind(th2));% Resulant velocity at y axis 
u=sqrt((ux.^2+uy.^2)); % That's Resultant Velocity at x and y 
P=P_max.*sqrt(1-((X./b).^2)-((X./b).^2));
gama=(wsp./w1);% Spin ratio 
J=((wsp.*sqrt(a.*b))./U1);% spin parameter
tau_epr_x=(32.*G./(pi.*(4-(3.*v))).*(2.*Cr./(2-Cr))).*((X+b)./(2.*b))% shear at pure rolling  
tau_epr_y=0;
tau_epr=sqrt(tau_epr_x.^2+tau_epr_y.^2);
tau_eps_x=(((64.*(2-v).*G)./(9.*pi.*(3-2.*v))).*J).*((Y./sqrt(a.*b)));
tau_eps_y=(((64.*(2-v).*G)./(9.*pi.*(3-2.*v))).*J).*((X./sqrt(a.*b)));
tau_esp=sqrt(tau_eps_x.^2+tau_eps_y.^2);
tau_e=abs(tau_epr+tau_esp);


%% so the calculate M (Traction cofficent)
figure(1)
I =trapz(y,trapz(x,tau_e,2))
Fc=linspace(0.1,2,100)
M=(abs((I./Fc)*10^-5));
% %% After calculations of the  values of Cr 
plot(linspace(0,0.02,100),-1*M)
title('Pure Rolling-Oil Santorac,Pmax=1.36 GPa,Temperature Effects Excluded')
xlabel('Creep Rate [C_r]')
ylabel(' Traction Coefficient[\mu]')
%% Plot Pressure 
figure(2)
plot(V_oil,P(1,:),'-r')
hold on 
plot(V_oil,P(2,:),'-r')
xlabel('Velocity of oil (m/s)')
ylabel('Pressure (Pa)')
%% clamping force
figure(3)
plot(P(1,:),Fc');
xlabel('Pressure (Pa)')
ylabel('Fc(clambing Force')
%% 
u1=70;    %linear velocity [m/s]
mumax=0.071;
Qmax=18800; %[kw] power
fc=Qmax/(mumax*u1);% clamp force calculated using Qmax. Please see page 5 under Test conditions
r1=70e-3;%[m]
j=0.0068;
gamma=0.29;
w1=u1/r1;%[rad/s]
wsp=w1*gamma;%spin rotation velocity
pmax=3e9;%[GPa]maximum pressure
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus
Rax=35e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=36.3e-3;%[m] roller 2
Rby=35.046e-3;% [m] roller 1
Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;
%coefficients to be found from book tribology by Batchelor
k0=sqrt((1/Rax-1/Ray)^2+(1/Rbx-1/Rby)^2+2*(1/Rax-1/Ray)*(1/Rbx-1/Rby)*cos(phi))/(1/Rax+1/Ray+1/Rbx+1/Rby);
k1=0.75; %values from graphs in tribology book by Batchelor
k2=0.61;
a=k1*(3*fc*Rprime/Eprime)^(1/3);
b=k2*(3*fc*Rprime/Eprime)^(1/3);
load 'ellipse-connectivity.txt'
ellipse = importdata('ellipse-connectivity.txt');
condata = ellipse(:,1);
load 'ellipse-nodes.txt'
Ellipse = importdata('ellipse-nodes.txt');
nodecoords=Ellipse(:,1);
[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation
%% 
Cr=0.04;
R1=140e-3%input('Rotational radius of drive side rolling element:');
R2=120e-3%input('Rotational radius of driven side rolling element:');
a=R1/8%input('Long radius of contact ellipse :');
b=R1/2%input('short radius of contact ellipse :');
U1=10%input('Peripheral Velocity of drive side rolling element : ');
w1=6%input('angular Velocity of drive side rolling element : ');
w2=3%input('angular Velocity of driven side rolling element : ');
U2=15%input('Peripheral Velocity of driven side rolling element : ');
Fc=20e3 %input('Clamping Force (Axial Force) : ' );
th1=45%input('Angle between tangent plane &Rotate axis of drive side rolling element :');
th2=36%input('Angle between tangent plane &Rotate axis of driven side rolling element :');
P_max =3e9%input('Max contact pressure  :');%   Max contact pressure
G=150e6%input(' Effictive Shear Modules :');
v=0.3%input('poison ratio of rolling element');
x = linspace(0,a/2,100); 
y = linspace(0,b/2,100); 
[X,Y] = meshgrid(x,y);
tau_e=[]
%% Calculates distributions of contact pressure & Velocity 
wsp= (w1.*sind(th1))-(w2.*sind(th2));% Spin Angular velocity 
ux=(U1.*cosd(th1))-(U2.*cosd(th2));% Resulant velocity  at x axis 
uy=(U1.*sind(th1))-(U2.*sind(th2));% Resulant velocity at y axis 
u=sqrt((ux.^2+uy.^2)); % That's Resultant Velocity at x and y 
P=P_max.*sqrt(1-((X./b).^2)-((X./b).^2));
gama=(wsp./w1);% Spin ratio 
J=((wsp.*sqrt(a.*b))./U1);% spin parameter
tau_epr_x=(32.*G./(pi.*(4-(3.*v))).*(2.*Cr./(2-Cr))).*((X+b)./(2.*b))% shear at pure rolling  
tau_epr_y=0;
tau_epr=sqrt(tau_epr_x.^2+tau_epr_y.^2);
tau_eps_x=(((64.*(2-v).*G)./(9.*pi.*(3-2.*v))).*J).*((Y./sqrt(a.*b)));
tau_eps_y=(((64.*(2-v).*G)./(9.*pi.*(3-2.*v))).*J).*((X./sqrt(a.*b)));
tau_esp=sqrt(tau_eps_x.^2+tau_eps_y.^2);
tau_e=tau_epr+tau_esp;
%% so the calculate M (Traction cofficent)
I = trapz(y,trapz(x,tau_e,2))
M=I/Fc;
clear all
 
%********************************************************************
%Choose roller dimensions and use working parameters to find clamping
%force Fc.
%********************************************************************
 
iN=60000;%input shaft RPM
iomega=iN*2*pi/60;%rad/s
ipower=15e3;%[W]
itorque=ipower/iomega;
oN=20000;%output shaft [RPM]
oomega=oN*2*pi/60;%rad/s
RR=iN/oN; %[reduction ratio]
 
%**********************************************************************
% design input shaft diameter based on torque and critical frequency. A
%factor of 5 was multiplied with shaftdia provide safety against critical 
%frequency and any other loading condition e.g. bending etc during application
%**********************************************************************
Su=2015;%%N/mm^2
ishaftdia=5*3*power((itorque*1000/Su),(1/3));%mm
YoungE=208e3;%N/mm^2
I=pi*ishaftdia^(4)/64;%mm^4
density=7.86e-3;%g/mm^3
A=pi*ishaftdia^(2)/4;%mm^2
g=9.81e3;%mm/s^2
l=50;%mm
iomegacritical=(pi/l)^(2)*sqrt((g*YoungE*I)/(A*density));
%************************************************************************
%use reduction ratio RR=2+2*(Dp/Ds) to get sun and planet diameters
%*************************************************************************
diasun=20; %[mm] choose as open choice
diaplanet=10;%[mm] use above RR equation
diaring=(diasun+2*diaplanet);%[mm]
otorque=RR*itorque;
 
% we have chosen 4 planet and sun configuration, therefore, using formula
% mud=Tr/(R*Fc)to find Fc, the clamping force
 
mud=0.075;%traction coefficient chosen 75% of the usual value for most lubricants
Fc=0.25*itorque/(mud*diasun*1e-3*0.5);%[N]
Nsun=iN;%[RPM]
Nplanet=Nsun/2;%[RPM]
radSunTrans=500;%[mm]
radPlanTrans=100;%[mm]
radRingTrans=100;
Rax=diasun/2*1e-3; %[m] 
Rbx=diaplanet/2*1e-3; % [m] 
Ray=radSunTrans*1e-3;%[m] 
Rby=radPlanTrans*1e-3;% [m] 
Rcx=-diaring/2*1e-3;
Rcy=radRingTrans*1e-3;
%**********************************************************************
%design output shaft
 
%**********************************************************************
oshaftdia=5*3*power((otorque*1000/Su),(1/3));%mm. A factor of 5 multiplied
% to make the shaft safe against whirling, to increase gap between
% rotational velocity and critical frequency, and to provide safety against
% unforeseen factors like material defects, machining related flaws
%and loading conditions
YoungE=208e3;%N/mm^2
I=pi*oshaftdia^(4)/64;%mm^4
density=7.86e-3;%g/mm^3
A=pi*oshaftdia^(2)/4;%mm^2
g=9.81e3;%mm/s^2
l=50;%mm
oomegacritical=(pi/l)^(2)*sqrt((g*YoungE*I)/(A*density));
 
%**********************************************************************
%Check for the required condition to choose the x and y direction radii
%can be activated to check
%**********************************************************************
% disp((1/Rbx+1/Rcx)>(1/Rby+1/Rcy));
% disp((1/Rax+1/Rbx)>(1/Ray+1/Rby));
 
%***********************************************************************
%Fatigue life calculations
%***********************************************************************
 
[rho]=findrfk(Rax,Rbx,Ray,Rby);
% 
% 
%                             % rho2=1/Rbx+1/Rcx+1/Rby+1/Rcy;
%                             % F2=(1/Rbx+1/Rcx-(1/Rby+1/Rcy))/rho2;
%                             % k22=1.175e6; %see figure 4
% 
k2=4.2e6;
Rsun=diasun/2*1e-3;
k4=2.32e19;%[N.m]
Lsun=k4*k2^0.9*Fc^(-3)*rho^(-6.3)*Rsun^(-0.9);%[millions of revolutions]
usun=4;%[number of sun stress cycles/rev]
SunLife=Lsun/(usun*Nsun)*1e6/60;%life in Hrs
Rplanet=diaplanet/2*1e-3;
Lplanet=k4*k2^0.9*Fc^(-3)*rho^(-6.3)*Rplanet^(-0.9);%[millions of revolutions]
uplanet=2;%[number of planet stress cycles/rev]
PlanetLife=Lplanet/(uplanet*Nplanet)*1e6/60;%life in Hrs
systemLife=((1/SunLife)^(10/9)+3*(1/PlanetLife)^(10/9))^(-0.9);%life in Hrs
 
                                % [rho,k2]=findrfk(Rbx,Rcx,Rby,Rcy);
                                % Rring=diaring/2*1e-3;
                                % Lring=k4*k2^0.9*Fc^(-3)*rho^(-6.3)*Rring^(-0.9);%[millions of revolutions]
                                % uring=4;
                                % RingLife=Lring/(uring*Nring)
 
%************************************************************************
%MAX DEFLECTION
%************************************************************************
nu=0.3;
E=208e9;
[delta1,a1,b1,pavg1]=finddelta(nu,E,Rax,Rbx,Ray,Rby,Fc);
[delta2,a2,b2,pavg2]=finddelta(nu,E,Rbx,Rcx,Rby,Rcy,Fc);
totaldeflection=2*(delta1+delta2);
disp('Total Deflection mm:');
fprintf('\n%f',totaldeflection*1e3);

%% 
% figure(4)
% temp = importdata('Temp-effect.txt');
% X = temp(:,1);
% Y = temp(:,2);
% hold on
% plot(X,Y,'o')
% x=0:0.01:0.5;
% y=0:0.0013:0.065;
% plot(x,y,'LineWidth',1.5)
% hold on
% x1=0.5:0.01:2;
% y1=0.058:0.0000465:0.065;
% yy=flipud(y1');
% plot(x1,yy,'LineWidth',1.5)
% xlim([0 2])
% ylim([0 0.1])
% title('Pure Rolling-Oil Santorac,U1=20.9 m/s,Pmax=1.36 GPa,Toil=140C, temperature effects included')
% xlabel('Creeping Rate in %age (Cr) ');
% ylabel('traction coefficient \mu');
% legend('Measured','Estimated')