function DCCON11(X,Y,D)
%       SUBROUTINE  DCCON1(X,Y,D)                                         09670000
%       IMPLICIT REAL*8 (A-H,O-Z)                                         09680000
% C                                                                       09690000
% C********************************************************************** 09700000
% C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** 09710000
% C********************************************************************** 09720000
% C                                                                       09730000
% C***** INPUT                                                            09740000
% C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    09750000
% C### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO  09760000
% C                                                                       09770000
%       COMMON /C0/DUMMY(5),SD,CD                                         09780000
%       COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,     09790000
%      *           UY,VY,WY,UZ,VZ,WZ                                      09800000
global DUMMY SD CD
global P Q S T XY X2 Y2 D2 R R2 R3 R5 QR QRX A3 A5 B3 C3 UY VY WY UZ VZ WZ
global N_CELL
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F3 = ones(N_CELL,1,'double').*3.0;
F5 = ones(N_CELL,1,'double').*5.0;
EPS = ones(N_CELL,1,'double').*0.000001;

%       DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,1.D-6/                  09810000
% C-----                                                                  09820000
c1 = abs(X) < EPS;   
c2 = abs(X) >= EPS;
X = F0.*c1 + X.*c2;
c1 = abs(Y) < EPS;
c2 = abs(Y) >= EPS;
Y = F0.*c1 + Y.*c2;
c1 = abs(D) < EPS;
c2 = abs(D) >= EPS;
D = F0.*c1 + D.*c2;
%       IF(DABS(X).LT.EPS) X=F0                                           09830000
%       IF(DABS(Y).LT.EPS) Y=F0                                           09840000
%       IF(DABS(D).LT.EPS) D=F0                                           09850000
	P=Y.*CD+D.*SD;  % p = y.*cos(dip) + d.*sin(dip)
	Q=Y.*SD-D.*CD;  % q = y.*sin(dip) - d.*cos(dip)
	S=P.*SD+Q.*CD;  % s = p.*sin(dip) + q.*cos(dip)
	T=P.*CD-Q.*SD;  % t = p.*cos(dip) - q.*sin(dip)
	XY=X.*Y;        %
	X2=X.*X;
	Y2=Y.*Y;
	D2=D.*D;
	R2=X2+Y2+D2;   % R^2=x^2 + y^2 + d^2
	R =sqrt(R2);
%       IF(R.EQ.F0) RETURN                                                09960000
c1 = R == F0;
if sum(rot90(sum(c1))) > 0
    return
end
    R3=R .*R2;
	R5=R3.*R2;
	R7=R5.*R2;
% C-----                                                                  10000000
	A3=F1-F3.*X2./R2;  % A3 = 1-3*x^2/R^2
	A5=F1-F5.*X2./R2;  % A5 = 1-5*x^2/R^2
	B3=F1-F3.*Y2./R2;  % B3 = 1-3*y^2/R^2
	C3=F1-F3.*D2./R2;  % C3 = 1-3*d^2/R^2
% C-----                                                                  10050000
	QR=F3.*Q./R5;      % qR = 3*q/R^5
	QRX=F5.*QR.*X./R2; % qRx = 5*qR*x/R^2
% C-----                                                                  10080000
	UY=SD-F5.*Y.*Q./R2;    % Uy = sin(dip) - 5*y*q/R^2
	UZ=CD+F5.*D.*Q./R2;    % Uz = cos(dip) + 5*d*q/R^2
	VY=S -F5.*Y.*P.*Q./R2; % Vy = s - 5*y*p*q/R^2
	VZ=T +F5.*D.*P.*Q./R2; % Vz = t + 5*d*p*q/R^2
	WY=UY+SD;              % Wy = Uy + sin(dip) = 2*sin(dip) - 5*y*q/R^2
	WZ=UZ+CD;              % Wz = Uz + cos(dip) = 2*cos(dip) + 5*d*q/R^2
%       RETURN                                                            10150000
%       END                                                               10160000
