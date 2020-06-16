function [Mo Mw] = cal_Mw(area,slip)
% This function is used for calculating the magnitude of a earthquake. Note
% that this function can only calculate the magnitude for only one sliding
% direction.
%
% ### Some basic knowledge about seismic moment(Mo) and moment magnitude(Mw)
%     Aki and Richards(2002) derived the following formula to calculate the
%     seismic moment(Mo) due to slip on the fault patch,
%         Mo  = u * L * W* D                        (1)
%where, u : the shear modulus.                   [Unit:Pa/Gpa]
%           In general, we often set u to 30GPa, or 3.0*10^10Pa.
%           1MPa = 10^6Pa. 1GPa = 10^9Pa. 1Bar = 10^5Pa
%           1Pa = 1*N/m^2 = 1*kg/(m.s^2)
%           1N = 10^5dyn
%       L : the length of fault patch.           [Unit:m]
%       W : the width of fault patch.            [Unit:m]
%       D : the dislocation/slip of fault patch. [Unit:m]
%       Mo: the seismic moment.                  [Unit:N.m]
%
%      Based on the seismic moment(Mo) we can get the moment magnitude(Mw),
%      Hanks and Kanamori(1979) derived the empirical formula to compute
%      the Mw based on Mo,
%          log10(10^7*Mo) = 3/2(Mw+10.7)            (2)
% namely,  Mw = 2*log10(Mo)/3 - 6.06                (3)
%
% Hanks and Kanamori.A moment magnitude scale.JGR(1979)
% Aki and Richard.Quantitative seismology, Sausalito, California:
%                 University Science Books(2002)
%
% Input 
%   area: the area of triangle patches.   [Unit:m^2]
%   slip: slip magnitude.                 [Unit:m]
% Output
%   Mo  : seismic moment.                 [Unit:N.m]
%   Mw  : moment magnitude.               [Unit:Mw]
% 
% By Shuai WANG @WHU 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% disp('    [Calculate Seismic Moment and Seismic Magnitude]');
km2m    = 1e3;
GPa2Pa  = 1e9;                                % GPa to Pa
shear_modulus = 33.*GPa2Pa;                   % shear modulus

area = area*km2m^2;
Mo = shear_modulus * sum(area.*slip);
Mw = 2/3 * log10(Mo) - 6.06;

format;
disp(['      ... the slip potency         is ',num2str(sum(area.*slip)),' m^3 ...']);
disp(['      ... the seismic moment(Mo)   is ',num2str(Mo),' N.m ...']);
disp(['      ... the moment magnitude(Mw) is ',num2str(Mw),' ...']);
end

