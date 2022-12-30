function [K11,K12,K21,K22,K11y,K12y,K21y,K22y,L1,L2] = GetDistMPCGains(Fb1,Gb11,Gb12,Qb1,Rb1,Fb2,Gb21,Gb22,Qb2,Rb2)
% compute MPC gains for two player game

Rt1 = Rb1 + Gb11'*Qb1*Gb11 + Gb21'*Qb2*Gb21;

% S1iy = Gbi1'*Qbi;
S11y = Gb11'*Qb1;
S12y = Gb21'*Qb2;


Rt2 = Rb2 + Gb22'*Qb2*Gb22 + Gb12'*Qb1*Gb12;

% S2iy = Gbi2'*Qbi;
S22y = Gb22'*Qb2;
S21y = Gb12'*Qb1;

% K1iy = Rt1^-1 * S1iy
K11y = Rt1^-1 * S11y;
K12y = Rt1^-1 * S12y;

% K2iy = Rt2^-1 * S2iy
K21y = Rt2^-1 * S21y;
K22y = Rt2^-1 * S22y;

%  K1i = -K1iy*Fbi
K11 = -K11y*Fb1; 
K12 =-K12y*Fb2; 

% K2i = -K2iy*Fbi
K21 = -K21y*Fb1;
K22 = -K22y*Fb2;


L1 = -K11y*Gb12 - K12y*Gb22;
L2 = -K22y*Gb21 - K21y*Gb11;
 