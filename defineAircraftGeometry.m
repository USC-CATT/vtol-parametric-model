function plane = defineAircraftGeometry()
%% Problem 3 calculating CoM and Inertia tensor for composite aircraft
%% Define Aircraft Mass and Geometry Properties (using grams)
% mass  xSize  ySize  zSize  xLoc   yLoc  zLoc
%    1  2      3      4      5      6     7
componentMassesAndGeom = ...
[   90  0.1    0.96   0.01  -0.23   0.44  0;    ... %  1 RightWing+Servo (s4)
    90  0.1    0.96   0.01  -0.23  -0.44  0;    ... %  2 LeftWing+Servo (s5)
%    13  0.075  0.35   0.002 -0.76   0    -0.16; ... % 3 Hor. Stab. (s2)OLD
   6.5  0.075  0.175  0.002 -0.76 -0.0875 -0.16; ... %  3 Hor. Stab. LFT(s2)
     0  0.08   0.002  0.18  -0.76   0    -0.09; ... %  4 Vert. Stab. (s3)
    72  0.065  0.035  0.015 -0.05   0     0.03; ... %  5 Battery
   106  0.87   0.07   0.07  -0.4    0     0;    ... %  6 Fuselage
    27  0.05   0.03   0.005 -0.05   0     0.02; ... %  7 Motor Controller
    10  0.04   0.02   0.005 -0.1    0     0.02; ... %  8 Radio
    20  0.05   0.01   0.01  -0.014  0     0;    ... %  9 2 Servos
    40  0.03   0.02   0.02   0.02   0     0.01; ... % 10 Motor
    12  0      0.26   0.025  0.05   0     0.01; ... % 11 Propeller
   6.5  0.075  0.175  0.002 -0.76 0.0875 -0.16];    % 12 Hor. Stab. RGT(s1)
n_masses = size(componentMassesAndGeom,1);

% Calculating Center of Mass
mass = 0;
x_CMrR_B = [0;0;0]; % location of center of mass relative to Reference location R
for i = 1:n_masses
    m = componentMassesAndGeom(i,1);
    dx = componentMassesAndGeom(i,5);
    dy = componentMassesAndGeom(i,6);
    dz = componentMassesAndGeom(i,7);
    
    mass = mass + m;
    x_CMrR_B = x_CMrR_B + m*[dx;dy;dz];

end

x_CMrR_B = x_CMrR_B/mass;
mass = mass/1000;

% Calculating Inertia Tensor
I_xx = 0;
I_yy = 0;
I_zz = 0;
I_xy = 0;
I_xz = 0;
I_yz = 0;

for i = 1:n_masses
    x = x_CMrR_B - componentMassesAndGeom(i,5:7)';
    m = componentMassesAndGeom(i,1);
    dx = componentMassesAndGeom(i,2:4)';
    I_xx = I_xx + 1/12*m*(dx(2)^2 + dx(3)^2) + m*(x(2)^2+x(3)^2);
    I_yy = I_yy + 1/12*m*(dx(1)^2 + dx(3)^2) + m*(x(1)^2+x(3)^2);
    I_zz = I_zz + 1/12*m*(dx(1)^2 + dx(2)^2) + m*(x(1)^2+x(2)^2);
    I_xy = I_xy + 0 + m*x(1)*x(2);
    I_xz = I_xz + 0 + m*x(1)*x(3);
    I_yz = I_yz + 0 + m*x(2)*x(3);
end

I_CM = [I_xx, I_xy, I_xz;...
        I_xy, I_yy, I_yz;...
        I_xz, I_yz, I_zz]/1000;
I_CM_inv = inv(I_CM);
plane.I_CM = I_CM;
plane.I_CM_inv = I_CM_inv;
plane.mass = mass;
plane.x_CMrR_B = x_CMrR_B;
fprintf('Total mass of vehicle is : %8.3f [kg]\n', mass)
fprintf('Center of mass relative to reference point is located at: \n')
fprintf('%8.3f [m]\n', x_CMrR_B')

fprintf('Inertia Tensor about center of mass is: \n')
fprintf('[%8.3f%8.3f%8.3f]\n[%8.3f%8.3f%8.3f] [kg-m^2]\n[%8.3f%8.3f%8.3f]\n', I_CM)
%% Problem 4, defining location of aero surfaces for force/moment calculation
n_s1 = [0; 0; -1]; %right side of horizontal stabilizer
n_s2 = [0; 0; -1]; 
n_s3 = [0; 1; 0]; 
n_s4 = [0; 0; -1]; 
n_s5 = [0; 0; -1];

CL0_s1 = 0; %right side of horizontal stabilizer
CL0_s2 = 0; 
CL0_s3 = 0; 
CL0_s4 = 0.05; 
CL0_s5 = 0.05;
e_s1 = 0.8;
e_s2 = 0.8; 
e_s3 = 0.8; 
e_s4 = 0.9; 
e_s5 = 0.9;
i_s1 = 0;
i_s2 = 0; 
i_s3 = 0; 
i_s4 = 0.05; 
i_s5 = 0.05;
CD0_s1 = 0.01;
CD0_s2 = 0.01; 
CD0_s3 = 0.01; 
CD0_s4 = 0.01; 
CD0_s5 = 0.01;
CDa_s1 = 1;
CDa_s2 = 1; 
CDa_s3 = 1; 
CDa_s4 = 1; 
CDa_s5 = 1;
a0_s1 = 0;
a0_s2 = 0; 
a0_s3 = 0; 
a0_s4 = 0.05; 
a0_s5 = 0.05;
CM0_s1 = 0;
CM0_s2 = 0; 
CM0_s3 = 0; 
CM0_s4 = -0.05; 
CM0_s5 = -0.05;
CMa_s1 = 0;
CMa_s2 = 0; 
CMa_s3 = 0; 
CMa_s4 = 0; 
CMa_s5 = 0;
% s1 Horizontal Stabilizer right side (Row 12 above)
c_s1 = componentMassesAndGeom( 12, 2);
b_s1 = componentMassesAndGeom( 12, 3);
S_s1 = b_s1*c_s1;
AR_s1 = b_s1/c_s1;
x_s1rR_B = [componentMassesAndGeom( 12, 5) + c_s1/4;... % position of s2 lift point relative to point R in body coordinates
            componentMassesAndGeom( 12, 6);...
            componentMassesAndGeom( 12, 7)];
% s2 Horizontal Stabilizer (Row 3 above)
c_s2 = componentMassesAndGeom( 3, 2);
b_s2 = componentMassesAndGeom( 3, 3);
S_s2 = b_s2*c_s2;
AR_s2 = b_s2/c_s2;
x_s2rR_B = [componentMassesAndGeom( 3, 5) + c_s2/4;... % position of s2 lift point relative to point R in body coordinates
            componentMassesAndGeom( 3, 6);...
            componentMassesAndGeom( 3, 7)];

% s3 Vertical Stabilizer (Row 4 above)
c_s3 = componentMassesAndGeom( 4, 2);
b_s3 = componentMassesAndGeom( 4, 4);   % this is the vertical stabilizer, so span is in the z dir.
S_s3 = b_s3*c_s3;
AR_s3 = b_s3/c_s3;
x_s3rR_B = [componentMassesAndGeom( 4, 5) + c_s3/4;... % position of s3 lift point relative to Reference location R in body coordinates
            componentMassesAndGeom( 4, 6);...
            componentMassesAndGeom( 4, 7)];

% s4 Right Wing (Row 1 above)
c_s4 = componentMassesAndGeom( 1, 2);
b_s4 = componentMassesAndGeom( 1, 3);
S_s4 = b_s4*c_s4;
AR_s4 = b_s4/c_s4;
x_s4rR_B = [componentMassesAndGeom( 1, 5) + c_s4/4;... % position of s4 lift point relative to Reference location R in body coordinates
            componentMassesAndGeom( 1, 6);...
            componentMassesAndGeom( 1, 7)];

% s5 Left Wing (Row 2 above)
c_s5 = componentMassesAndGeom( 2, 2);
b_s5 = componentMassesAndGeom( 2, 3);
S_s5 = b_s5*c_s5;
AR_s5 = b_s5/c_s5;
x_s5rR_B = [componentMassesAndGeom( 2, 5) + c_s5/4;... % position of s5 lift point relative to Reference location R in body coordinates
            componentMassesAndGeom( 2, 6);...
            componentMassesAndGeom( 2, 7)];


fprintf(' =======================================================================\n')
fprintf('|         Aerodynamic Locations and Properties of Aero Surfaces         |\n')
fprintf('|________|_________x_s#/R^B_________|___________________________________|\n')
fprintf('|___s#___|____x___|____y___|____z___|____c___|____b___|____S___|___AR___|\n')
fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
    's1',x_s1rR_B, c_s1, b_s1, S_s1, AR_s1)
fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
    's2',x_s2rR_B, c_s2, b_s2, S_s2, AR_s2)
fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
    's3',x_s3rR_B, c_s3, b_s3, S_s3, AR_s3)
fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
    's4',x_s4rR_B, c_s4, b_s4, S_s4, AR_s4)
fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
    's5',x_s5rR_B, c_s5, b_s5, S_s5, AR_s5)
fprintf(' =======================================================================\n')

%% Problem 5, compiling aerodynamic terms for use in runtime.
% position of aero surfaces WRT CM
plane.x_srCM_B = [x_s1rR_B - x_CMrR_B,...
            x_s2rR_B - x_CMrR_B,...
            x_s3rR_B - x_CMrR_B,...
            x_s4rR_B - x_CMrR_B,...
            x_s5rR_B - x_CMrR_B];


plane.x_RrCM_B = [plane.x_srCM_B(:,4:5), [componentMassesAndGeom(6,5);...
    componentMassesAndGeom(6,6); componentMassesAndGeom(6,7)]-x_CMrR_B]
%% combine all variables into single variables for ease of iteration
plane.n   = [  n_s1,   n_s2,   n_s3,   n_s4,   n_s5];
plane.CL0 = [CL0_s1, CL0_s2, CL0_s3, CL0_s4, CL0_s5];
plane.e   = [  e_s1,   e_s2,   e_s3,   e_s4,   e_s5];
plane.i_s = [  i_s1,   i_s2,   i_s3,   i_s4,   i_s5];
plane.CD0 = [CD0_s1, CD0_s2, CD0_s3, CD0_s4, CD0_s5];
plane.CDa = [CDa_s1, CDa_s2, CDa_s3, CDa_s4, CDa_s5];
plane.a0  = [ a0_s1,  a0_s2,  a0_s3,  a0_s4,  a0_s5];
plane.CM0 = [CM0_s1, CM0_s2, CM0_s3, CM0_s4, CM0_s5];
plane.CMa = [CMa_s1, CMa_s2, CMa_s3, CMa_s4, CMa_s5];
plane.c   = [  c_s1,   c_s2,   c_s3,   c_s4,   c_s5];
plane.b   = [  b_s1,   b_s2,   b_s3,   b_s4,   b_s5];
plane.S   = [  S_s1,   S_s2,   S_s3,   S_s4,   S_s5];
plane.AR  = [ AR_s1,  AR_s2,  AR_s3,  AR_s4,  AR_s5];
plane.n_aero = 5; % number of aero surfaces

plane.rho = 1.225;    % [kg/m^3] density of air at sea level

plane.halfRhoS = 1/2*plane.rho*plane.S;
plane.halfRhoSc = plane.halfRhoS.*plane.c;
plane.CLa = 2*pi*plane.AR./(2+plane.AR);
end
