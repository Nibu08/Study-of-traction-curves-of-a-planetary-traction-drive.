
%% OIl III 2.36 %*************************************************************
%Given geometry and material data
%*************************************************************
clc
clear all
u1=20.9;    %linear velocity [m/s]
a=0.001; % assumed contact elliptical area dimensions,
%later corrected thru iteration. see line 49 of this code
b=0.001;
Qmax=21000; %[kw] power
j=0;
r1=70e-3;%[m]
gamma=0;
w1=u1/r1;%[rad/s]
wsp=w1*gamma;%spin rotation velocity
pmax=2.36e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus
 
%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************
 
%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1
 
Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact
 
%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;
 
% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
fc=2*pi*a*b/3*pmax;
a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
tmpa=a;
tmpb=b;
disp(2*a);
disp(2*b);
disp(fc);
end
 
%**********************************************************
%preprocessing for contact point area geometry
 
%**********************************************************
fid=fopen('ellipse-connectivity-oil-A-III.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);
 
fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);
 
[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation
 
for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);
        
    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;
 
xx=centroidx;
yy=centroidy;
 
%**********************************************************************
 
%calculate shear stress distribution
 
%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas
 
k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;
 
for cr=0:0.001:0.02
 
teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
tepry=0;
 
tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);
 
te=teprx;
 
for i=1:q
    
    if te(i)>=tl(i)
        te(i)=tl(i);
    end
end
 
ftraction=sum(te.*area);
 
mu(k)=ftraction/fc;
k=k+1;
end
cr=0:0.001:0.02;
figure (6)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title('Pure Rolling, Oil A-III, U1=20.9 m/s,Pmax=3.5 Gpa,T=80')
xlabel('Creeping Rate in %age (Cr) ');
ylabel('traction coefficient \mu');
ylim([0 0.08])

%% %% %*************************************************************
%Given geometry and material data
%*************************************************************
clc
clear all
u1=20.9;    %linear velocity [m/s]
a=0.001; % assumed contact elliptical area dimensions,
%later corrected thru iteration. see line 49 of this code
b=0.001;
Qmax=21000; %[kw] power
j=0;
r1=70e-3;%[m]
gamma=0;
w1=u1/r1;%[rad/s]
wsp=w1*gamma;%spin rotation velocity
pmax=3.18e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus
 
%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************
 
%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1
 
Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact
 
%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;
 
% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
fc=2*pi*a*b/3*pmax;
a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
tmpa=a;
tmpb=b;
disp(2*a);
disp(2*b);
disp(fc);
end
 
%**********************************************************
%preprocessing for contact point area geometry
 
%**********************************************************
fid=fopen('ellipse-connectivity-oil-A-III.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);
 
fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);
 
[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation
 
for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);
        
    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;
 
xx=centroidx;
yy=centroidy;
 
%**********************************************************************
 
%calculate shear stress distribution
 
%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas
 
k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;
 
for cr=0:0.001:0.02
 
teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
tepry=0;
 
tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);
 
te=teprx;
 
for i=1:q
    
    if te(i)>=tl(i)
        te(i)=tl(i);
    end
end
 
ftraction=sum(te.*area);
 
mu(k)=ftraction/fc;
k=k+1;
end
cr=0:0.001:0.02;
figure (7)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title('Pure Rolling, Oil A-III, U1=20.9 m/s,Pmax=3.18 Gpa')
xlabel('Creeping Rate in %age (Cr) ');
ylabel('traction coefficient \mu');
xlim([0 0.014])
ylim([0 0.04])
%% %% %% %*************************************************************
%Given geometry and material data
%*************************************************************
clc
clear all
u1=30;    %linear velocity [m/s]
a=0.001; % assumed contact elliptical area dimensions,
%later corrected thru iteration. see line 49 of this code
b=0.001;
Qmax=21000; %[kw] power
j=0;
r1=70e-3;%[m]
gamma=0;
w1=u1/r1;%[rad/s]
wsp=w1*gamma;%spin rotation velocity
pmax=2.8e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus
 
%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************
 
%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1
 
Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact
 
%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;
 
% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
fc=2*pi*a*b/3*pmax;
a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
tmpa=a;
tmpb=b;
disp(2*a);
disp(2*b);
disp(fc);
end
 
%**********************************************************
%preprocessing for contact point area geometry
 
%**********************************************************
fid=fopen('ellipse-connectivity-oil-A-III.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);
 
fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);
 
[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation
 
for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);
        
    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;
 
xx=centroidx;
yy=centroidy;
 
%**********************************************************************
 
%calculate shear stress distribution
 
%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas
 
k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;
 
for cr=0:0.001:0.02
 
teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
tepry=0;
 
tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);
 
te=teprx;
 
for i=1:q
    
    if te(i)>=tl(i)
        te(i)=tl(i);
    end
end
 
ftraction=sum(te.*area);
 
mu(k)=ftraction/fc;
k=k+1;
end
cr=0:0.001:0.02;
figure (8)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title('Pure Rolling, Oil A-III, U1=10 m/s,Pmax=2.3 Gpa,T=140')
xlabel('Creeping Rate in %age (Cr) ');
ylabel('traction coefficient \mu');
xlim([0 0.014])
ylim([0 0.06])
%% 
clc
clear all
u1=20;    %linear velocity [m/s]
a=0.001; % assumed contact elliptical area dimensions,
%later corrected thru iteration. see line 49 of this code
b=0.001;
Qmax=21000; %[kw] power
j=0;
r1=70e-3;%[m]
gamma=0;
w1=u1/r1;%[rad/s]
wsp=w1*gamma;%spin rotation velocity
pmax=2.5e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus
 
%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************
 
%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1
 
Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact
 
%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;
 
% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
fc=2*pi*a*b/3*pmax;
a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
tmpa=a;
tmpb=b;
disp(2*a);
disp(2*b);
disp(fc);
end
 
%**********************************************************
%preprocessing for contact point area geometry
 
%**********************************************************
fid=fopen('ellipse-connectivity-oil-A-III.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);
 
fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);
 
[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation
 
for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);
        
    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;
 
xx=centroidx;
yy=centroidy;
 
%**********************************************************************
 
%calculate shear stress distribution
 
%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas
 
k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;
 
for cr=0:0.001:0.02
 
teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
tepry=0;
 
tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);
 
te=teprx;
 
for i=1:q
    
    if te(i)>=tl(i)
        te(i)=tl(i);
    end
end
 
ftraction=sum(te.*area);
 
mu(k)=ftraction/fc;
k=k+1;
end
cr=0:0.001:0.02;
figure (9)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title('Pure Rolling, Oil A-III, U1=20 m/s,Pmax=2.5 Gpa,Toil=120C, Temperature Effects Included')
xlabel('Creeping Rate in %age (Cr) ');
ylabel('traction coefficient \mu');
ylim([0 0.06])


