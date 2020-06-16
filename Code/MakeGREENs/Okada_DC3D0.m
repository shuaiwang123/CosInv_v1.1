% function[UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET] = Okada_DC3D0(ALPHA,...
%                 X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4)
function [UX,UY,UZ] = Okada_DC3D0(ALPHA,...
                X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4)
            
% 	SUBROUTINE  DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,      00010000
%      *               UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET) 00020002
% 	IMPLICIT REAL*8 (A-H,O-Z)                                         00030000
% 	REAL*4   ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,               00040000
%      *         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             00050000
% C                                                                       00060000
% C********************************************************************   00070000
% C*****                                                          *****   00080000
% C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   00090000
% C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   00100000
% C*****                         CODED BY  Y.OKADA ... SEP.1991   *****   00110002
% C*****                         REVISED   Y.OKADA ... NOV.1991   *****   00120002
% C*****                                                          *****   00130000
% C********************************************************************   00140000
% C                                                                       00150000
% C***** INPUT                                                            00160000
% C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           00170000
% C*****   X,Y,Z : COORDINATE OF OBSERVING POINT in the fault system      00180000
% C*****   DEPTH : SOURCE DEPTH                                           00190000
%                  the unit of X,Y,Z and Depth are given in the unit of m
% C*****   DIP   : DIP-ANGLE (DEGREE)                                     00200000
% C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        00210000
% C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2   00220000
% C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3     00230000
% C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4     00240000
%              the unit of Potency is given in the unit of m^3
% C                                                                       00250000
% C***** OUTPUT                                                           00260000
% C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /          00270000
% C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )  00280000
% C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /          00290000
% C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )  00300000
% C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     00310000
% C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  00320002
% C                                                                       00330000

    global DUMMY R
	global N_CELL
% number of observation points
    N_CELL = numel(X);
    
% 	COMMON /C1/DUMMY(8),R                                             00340000
% 	DIMENSION  U(12),DUA(12),DUB(12),DUC(12)                          00350000
% 	DATA  F0/0.D0/   

F0 = zeros(N_CELL,1,'double');
if Z>0.0
    disp(' ** POSITIVE Z WAS GIVEN IN SUB-DC3D');
end
        U  (1:N_CELL,1:12)=0.0;
        DUA(1:N_CELL,1:12)=0.0;
        DUB(1:N_CELL,1:12)=0.0;
        DUC(1:N_CELL,1:12)=0.0;
AALPHA=ALPHA;
DDIP=DIP;

DCCON0(AALPHA,DDIP);

% C======================================                                 00480000
% C=====  REAL-SOURCE CONTRIBUTION  =====                                 00490000
% C======================================                                 00500000
	XX=X;
	YY=Y;
	ZZ=Z;
	DD=DEPTH+Z;

DCCON1(XX,YY,DD);

% C=======================================                                00960000
% C=====  IN CASE OF SINGULAR (R=0)  =====                                00970000
% C=======================================
c1 = R == F0;
if sum(rot90(sum(c1))) > 1
        UX=F0;
        UY=F0;
        UZ=F0;
        UXX=F0;
        UYX=F0;
        UZX=F0;
        UXY=F0;
        UYY=F0;
        UZY=F0;
        UXZ=F0;
        UYZ=F0;
        UZZ=F0;
        IRET=ones(N_CELL,1);           
    return
end

	PP1=POT1; % Potency for strike-slip fault =(moment of double-couple)/u =LWU
	PP2=POT2; % Potency for dip-slip fault =(moment of double-couple)/u =LWU
	PP3=POT3; % Potency for tensile fault = (moment of isotropic part of dipole)/lamada
	PP4=POT4; % Potency for explosive source = (moment of dipole)/u

DUA = UA0(XX,YY,DD,PP1,PP2,PP3,PP4);
for I = 1:12
    if(I<10) 
        U(:,I)=U(:,I)-DUA(:,I);
    end
    if(I>=10)
        U(:,I)=U(:,I)+DUA(:,I);
    end
end

% C=======================================                                00670000
% C=====  IMAGE-SOURCE CONTRIBUTION  =====                                00680000
% C=======================================                                00690000
	DD=DEPTH-Z;
    DCCON1(XX,YY,DD);
    DUA = UA0(XX,YY,DD,PP1,PP2,PP3,PP4);
    DUB = UB0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4);
    DUC = UC0(XX,YY,DD,ZZ,PP1,PP2,PP3,PP4);
    for I = 1:12
        DU = DUA(:,I)+DUB(:,I)+ZZ.*DUC(:,I);
        if I>=10
            DU=DU+DUC(:,I-9);
        end
        U(:,I)=U(:,I)+DU;
    end
% C=====                                                                  00810000
% ./(1.0E+6) ---> unit adjustment. Strain is going to be further divided
%  by e+3 later
%   Area: the unit of the point source area is given in the unit of m^2
%   Slip: the unit of the slip is given in the unit of m
%
%         UX=U(:,1)./(1.0E+6);
%         UY=U(:,2)./(1.0E+6);
%         UZ=U(:,3)./(1.0E+6);
%         
%         UXX=U(:,4)./(1.0E+9);
%         UYX=U(:,5)./(1.0E+9);
%         UZX=U(:,6)./(1.0E+9);
%         
%         UXY=U(:,7)./(1.0E+9);
%         UYY=U(:,8)./(1.0E+9);
%         UZY=U(:,9)./(1.0E+9);
%         
%         UXZ=U(:,10)./(1.0E+9);
%         UYZ=U(:,11)./(1.0E+9);
%         UZZ=U(:,12)./(1.0E+6);
        UX=U(:,1)./(1.0E+0);
        UY=U(:,2)./(1.0E+0);
        UZ=U(:,3)./(1.0E+0);
        
        UXX=U(:,4)./(1.0E+0);
        UYX=U(:,5)./(1.0E+0);
        UZX=U(:,6)./(1.0E+0);
        
        UXY=U(:,7)./(1.0E+0);
        UYY=U(:,8)./(1.0E+0);
        UZY=U(:,9)./(1.0E+9);
        
        UXZ=U(:,10)./(1.0E+0);
        UYZ=U(:,11)./(1.0E+0);
        UZZ=U(:,12)./(1.0E+0);
        IRET=F0;
       
%         % [define Greens function]
%         % define the number of surface displacement directions(3=ENU)
%         n_surface_displacement_directions = 3;
%         % define the number of slip displacement directions (2=SS,DS)
%         n_slip_directions = 2;
%         %
        
        
        
        
        
        
        
