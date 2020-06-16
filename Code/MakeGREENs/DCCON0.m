function DCCON0(ALPHA,DIP)
% Okada 92 code subroutine DCCON0
%
global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D
global N_CELL
%       DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/             %09430000
%       DATA EPS/1.D-6/ 
F0 = zeros(N_CELL,1,'double');
F1 = ones(N_CELL,1,'double');
F2 = ones(N_CELL,1,'double').*2.0;
PI2 = ones(N_CELL,1,'double').*6.283185307179586;
EPS = ones(N_CELL,1,'double').*1.0e-6;

      ALP1=(F1-ALPHA)./F2;    % (1-alpha)./2
      ALP2= ALPHA./F2;        % alpha./2
      ALP3=(F1-ALPHA)./ALPHA; % (1-alpha)./alpha
      ALP4= F1-ALPHA;         % 1-alpha
      ALP5= ALPHA;            % alpha

      P18=PI2./double(360.0);                                                    %09520000
      SD=double(sin(DIP.*P18));                                                  %09530000
      CD=double(cos(DIP.*P18)); 
        c1 = abs(CD) < EPS;     % cos(dip)¡Ö0   dip¡Ö90, c1
        c2 = abs(CD) >= EPS;    % cos(dip)>0           , c2=1
        s1 = SD > F0;           % sin(dip)>0 
        s2 = SD == F0;          % sin(dip)=0 dip=0
        s3 = SD < F0;           % sin(dip)<0 

        CD = F0.*c1 + CD.*c2;   % cos(dip)=

% CAUTION ************ modified by Shinji Toda (CD = 0.0 produces 'nan')
%                      in MATLAB
%         c3 = abs(CD) < EPS;
%         c4 = abs(CD) <= EPS;
%         CD = c3.*EPS + c4.*CD;
% CAUTION ***************************************************************

%09560000
%     if SD>F0
%        SD= F1;
%     end
%     if SD<F0
%         SD=-F1;                                             %09580000
%     end
    SD = c1.*(F1.*s1 + SD.*s2 + (-1.0).*F1.*s3) + c2.*SD;% 
%end
                                                            %09590000
      SDSD=SD.*SD;            % sin(dip)^2                                           %09600000
      CDCD=CD.*CD;            % cos(dip)^2                                           %09610000
      SDCD=SD.*CD;            % sin(dip).*cos(dip)                                   %09620000
      S2D=F2.*SDCD;           % sin(2*dip)=2*sin(dip).*cos(dip)                      %09630000
      C2D=CDCD-SDSD;          % cos(2*dip)=cos(dip)*cos(dip)-sin(dip)*sin(dip)       %09640000
%       RETURN                                                            %09650000
%       END                                                               %09660000
