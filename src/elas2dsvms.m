function [Sed,Vms]=elas2dsvms(p,t,uu,E,nu) %-----------------------------------------
% Two-dimensional linear elasticity
% Sed Shear energy density
% Vms Von Mises effective stress
%----------------------------------------
n=size(p,2); nn=2*n;
% Lame constant
lambda=E*nu/((1+nu)*(1-2*nu)); mu=E/(2*(1+nu));
% area of traingles
it1=t(1,:); it2=t(2,:); it3=t(3,:);
x21=(p(1,it2)-p(1,it1))'; x32=(p(1,it3)-p(1,it2))'; x31=(p(1,it3)-p(1,it1))'; 
y21=(p(2,it2)-p(2,it1))'; y32=(p(2,it3)-p(2,it2))'; y31=(p(2,it3)-p(2,it1))'; 
ar=.5*(x21.*y31-x31.*y21);
% gradient of scalar basis functions
phi1=[-y32./(2*ar) x32./(2*ar)];
phi2=[ y31./(2*ar) -x31./(2*ar)];
phi3=[-y21./(2*ar) x21./(2*ar)];
% displacements
u=uu(1:2:end); v=uu(2:2:end);
uh=[u(it1) u(it2) u(it3)]; vh=[v(it1) v(it2) v(it3)];
% strains
e11=uh(:,1).*phi1(:,1)+uh(:,2).*phi2(:,1)+uh(:,3).*phi3(:,1); 
e22=vh(:,1).*phi1(:,2)+vh(:,2).*phi2(:,2)+vh(:,3).*phi3(:,2); 
e12=uh(:,1).*phi1(:,2)+uh(:,2).*phi2(:,2)+uh(:,3).*phi3(:,2)...
	+vh(:,1).*phi1(:,1)+vh(:,2).*phi2(:,1)+vh(:,3).*phi3(:,1); 
clear uh vh
% stresses
sig11=(lambda+2*mu)*e11+lambda*e22; sig22=lambda*e11+(lambda+2*mu)*e22; sig12=mu*e12;
clear e11 e22 e12
% area of patches 
arp=full(sparse(it1,1,ar,n,1)+sparse(it2,1,ar,n,1)+sparse(it3,1,ar,n,1));
% mean value of stresses on patches
sm1=ar.*sig11; sm2=ar.*sig22; sm12=ar.*sig12; 
s1=full(sparse(it1,1,sm1,n,1)+sparse(it2,1,sm1,n,1)+sparse(it3,1,sm1,n,1)); 
s2=full(sparse(it1,1,sm2,n,1)+sparse(it2,1,sm2,n,1)+sparse(it3,1,sm2,n,1)); 
s12=full(sparse(it1,1,sm12,n,1)+sparse(it2,1,sm12,n,1)+sparse(it3,1,sm12,n,1)); 
s1=s1./arp; s2=s2./arp; s12=s12./arp;
% Shear energy density 
Sed=((.5+mu*mu/(6*(mu+lambda)^2))*(s1+s2).^2+2*(s12.^2-s1.*s2))/(4*mu);
% Von Mises effective stress
delta=sqrt((s1-s2).^2+4*s12.^2); sp1=s1+s2+delta; sp2=s1+s2-delta; 
Vms=sqrt(sp1.^2+sp2.^2-sp1.*sp2);