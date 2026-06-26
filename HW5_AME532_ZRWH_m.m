%% HW5
clear variables

%% Preamble
t0 = 0;
tf = 2;
dt = 0.1;

%% Problem 1, gusty wind instantaneous angle of attack
v_BrW_W = [500;0;0]; % [ft/s] airspeed
a = deg2rad(8); % [rad] starting alpha
b = deg2rad(-5); % [rad] starting beta

C_WfB = [ cos(a)*cos(b)  sin(b)  sin(a)*cos(b);...
         -cos(a)*sin(b)  cos(b) -sin(a)*sin(b);...
            -sin(a)        0        cos(a)];
C_BfW = C_WfB';
v_BrW_B = C_BfW*v_BrW_W;
fprintf('Problem 1:\n')

% part a
v_GrB_B = [0;20;0]; % [ft/s] wind gust
v_BrWpG_B = v_BrW_B - v_GrB_B;
V_T = sqrt(v_BrWpG_B'*v_BrWpG_B);
b_a = asin(v_BrWpG_B(2)/V_T);
a_a = asin(v_BrWpG_B(3)/V_T/cos(b_a));
fprintf('Part a\n')
fprintf('\tStarting Airspeed in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_BrW_B)
fprintf('\talpha = %7.2f  [deg]\n\tbeta  = %7.2f  [deg]\n',...
    rad2deg(a),rad2deg(b))
fprintf('\tWind Gust in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_GrB_B)
fprintf('\tResultant Airspeed in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_BrWpG_B)
fprintf('    NEW alpha = %7.2f  [deg]\n    NEW beta  = %7.2f  [deg]\n',...
    rad2deg(a_a),rad2deg(b_a))
% part b
v_GrB_B = [50;0;0]; % [ft/s] wind gust
v_BrWpG_B = v_BrW_B - v_GrB_B;
V_T = sqrt(v_BrWpG_B'*v_BrWpG_B);
b_b = asin(v_BrWpG_B(2)/V_T);
a_b = asin(v_BrWpG_B(3)/V_T/cos(b_b));

fprintf('Part b\n')
fprintf('\tStarting Airspeed in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_BrW_B)
fprintf('\talpha = %7.2f  [deg]\n\tbeta  = %7.2f  [deg]\n',...
    rad2deg(a),rad2deg(b))
fprintf('\tWind Gust in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_GrB_B)
fprintf('\tResultant Airspeed in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_BrWpG_B)
fprintf('    NEW alpha = %7.2f  [deg]\n    NEW beta  = %7.2f  [deg]\n',...
    rad2deg(a_b),rad2deg(b_b))
% part c
v_GrB_B = [0;-30*cos(deg2rad(70));-30*sin(deg2rad(70))]; % [ft/s] wind gust
% not sure if you meant the gust was coming from below or going below the
% xy plane...
v_BrWpG_B = v_BrW_B - v_GrB_B;
V_T = sqrt(v_BrWpG_B'*v_BrWpG_B);
b_c = asin(v_BrWpG_B(2)/V_T);
a_c = asin(v_BrWpG_B(3)/V_T/cos(b_c));
fprintf('Part c\n')
fprintf('\tStarting Airspeed in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_BrW_B)
fprintf('\talpha = %7.2f  [deg]\n\tbeta  = %7.2f  [deg]\n',...
    rad2deg(a),rad2deg(b))
fprintf('\tWind Gust in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_GrB_B)
fprintf('\tResultant Airspeed in body frame:\n')
fprintf('\t\t%8.3f [ft/s]\n',v_BrWpG_B)
fprintf('    NEW alpha = %7.2f  [deg]\n    NEW beta  = %7.2f  [deg]\n',...
    rad2deg(a_c),rad2deg(b_c))

%% Problem 2, verifying relative wind calculation
fprintf('Case 1: Traveling straight forward, no side slip, no beta, no wind\n')
v_BrE_B = [100;0;0];
v_WrE_N = [0;0;0];
Phi_BrN = [0;0;0];
sim("HW5p2_AME532_ZRWH_slx")

fprintf('Inputs:\n')
fprintf('\tBody velocity:\n')
fprintf('\t\t%5.2f\n',v_BrE_B)
fprintf('\tWind velocity:\n')
fprintf('\t\t%5.2f\n',v_WrE_N)
fprintf('\tBody attitude:\n')
fprintf('\t\t%5.2f\n',Phi_BrN)
fprintf(['Expected results:\n',...
'\tv_BrW_B = [100;0;0];\n',...
'\talpha = 0;\n',...
'\tbeta = 0;\n'])
fprintf('Actual Results:\n')
fprintf('\tv_BrW_B = ;\n')
fprintf('\t\t%5.2f\n',v_BrW_B(1,:))
fprintf('\talpha = %5.2f;\n',alpha(1))
fprintf('\tbeta = %5.2f;\n',beta(1))

%}
%%% 
fprintf('Case 2: Traveling 30 degrees side slip, no wind\n')
v_BrE_B = [100*cos(pi/6);100*sin(pi/6);0];
v_WrE_N = [0;0;0];
Phi_BrN = [0;0;0];
sim("HW5p2_AME532_ZRWH_slx")

fprintf('Inputs:\n')
fprintf('\tBody velocity:\n')
fprintf('\t\t%5.2f\n',v_BrE_B)
fprintf('\tWind velocity:\n')
fprintf('\t\t%5.2f\n',v_WrE_N)
fprintf('\tBody attitude:\n')
fprintf('\t\t%5.2f\n',Phi_BrN)
fprintf(['Expected results:\n',...
'\tv_BrW_B = [86.6;50;0];\n',...
'\talpha = 0;\n',...
'\tbeta = 0.5236;\n'])
%{
Expected results:
v_BrW_B = [86.6;50;0];
alpha = 0;
beta = 0.5236; % pi/6
%}
fprintf('Actual Results:\n')
fprintf('\tv_BrW_B = ;\n')
fprintf('\t\t%5.2f\n',v_BrW_B(1,:))
fprintf('\talpha = %5.2f;\n',alpha(1))
fprintf('\tbeta = %5.2f;\n',beta(1))

%%% 
fprintf('Case 3: Traveling 30 degrees down, no wind\n')
v_BrE_B = [100*cos(pi/6);0;100*sin(pi/6)];
v_WrE_N = [0;0;0];
Phi_BrN = [0;0;0];
sim("HW5p2_AME532_ZRWH_slx")

fprintf('Inputs:\n')
fprintf('\tBody velocity:\n')
fprintf('\t\t%5.2f\n',v_BrE_B)
fprintf('\tWind velocity:\n')
fprintf('\t\t%5.2f\n',v_WrE_N)
fprintf('\tBody attitude:\n')
fprintf('\t\t%5.2f\n',Phi_BrN)
fprintf(['Expected results:\n',...
'\tv_BrW_B = [86.6;0;50];\n',...
'\talpha = 0.52;\n',...
'\tbeta = 0;\n'])
%{
Expected results:
v_BrW_B = [86.6;0;50];
alpha = 0.5236; % pi/6
beta = 0; %
%}
fprintf('Actual Results:\n')
fprintf('\tv_BrW_B = ;\n')
fprintf('\t\t%5.2f\n',v_BrW_B(1,:))
fprintf('\talpha = %5.2f;\n',alpha(1))
fprintf('\tbeta = %5.2f;\n',beta(1))
%%% 
fprintf('Case 4: Traveling straight forward, cross wind\n')
v_BrE_B = [86;0;0];
v_WrE_N = [0;50;0];
Phi_BrN = [0;0;0];
sim("HW5p2_AME532_ZRWH_slx")

fprintf('Inputs:\n')
fprintf('\tBody velocity:\n')
fprintf('\t\t%5.2f\n',v_BrE_B)
fprintf('\tWind velocity:\n')
fprintf('\t\t%5.2f\n',v_WrE_N)
fprintf('\tBody attitude:\n')
fprintf('\t\t%5.2f\n',Phi_BrN)
fprintf(['Expected results:\n',...
'\tv_BrW_B = [86;-50;0];\n',...
'\talpha = 0;\n',...
'\tbeta = -0.52;\n'])
%{
Expected results:
v_BrW_B = [86;-50;0];
alpha = 0; % 
beta = -0.5236; % pi/6
%}
fprintf('Actual Results:\n')
fprintf('\tv_BrW_B = ;\n')
fprintf('\t\t%5.2f\n',v_BrW_B(1,:))
fprintf('\talpha = %5.2f;\n',alpha(1))
fprintf('\tbeta = %5.2f;\n',beta(1))
%% Problems 3-5
% see defineAircraftGeometry.m
plane = defineAircraftGeometry();
busInfo = Simulink.Bus.createObject(plane);
% % Define Aircraft Mass and Geometry Properties (using grams)
% mass  xSize  ySize  zSize  xLoc   yLoc  zLoc
%    1  2      3      4      5      6     7
% componentMassesAndGeom = ...
% [   90  0.1    0.96   0.01  -0.23   0.44  0;    ... %  1 RightWing+Servo (s4)
%     90  0.1    0.96   0.01  -0.23  -0.44  0;    ... %  2 LeftWing+Servo (s5)
%     13  0.075  0.35   0.002 -0.76   0    -0.16; ... %  3 Hor. Stab. (s2)
%      0  0.08   0.002  0.18  -0.76   0    -0.09; ... %  4 Vert. Stab. (s3)
%     72  0.065  0.035  0.015 -0.05   0     0.03; ... %  5 Battery
%    106  0.87   0.07   0.07  -0.4    0     0;    ... %  6 Fuselage
%     27  0.05   0.03   0.005 -0.05   0     0.02; ... %  7 Motor Controller
%     10  0.04   0.02   0.005 -0.1    0     0.02; ... %  8 Radio
%     20  0.05   0.01   0.01  -0.014  0     0;    ... %  9 2 Servos
%     40  0.03   0.02   0.02   0.02   0     0.01; ... % 10 Motor
%     12  0      0.26   0.025  0.05   0     0.01];    % 11 Propeller
% n_masses = size(componentMassesAndGeom,1);
% 
% Calculating Center of Mass
% mass = 0;
% x_CMrR_B = [0;0;0]; % location of center of mass relative to Reference location R
% for i = 1:n_masses
%     m = componentMassesAndGeom(i,1);
%     dx = componentMassesAndGeom(i,5);
%     dy = componentMassesAndGeom(i,6);
%     dz = componentMassesAndGeom(i,7);
% 
%     mass = mass + m;
%     x_CMrR_B = x_CMrR_B + m*[dx;dy;dz];
% 
% end
% 
% x_CMrR_B = x_CMrR_B/mass;
% 
% Calculating Inertia Tensor
% I_xx = 0;
% I_yy = 0;
% I_zz = 0;
% I_xy = 0;
% I_xz = 0;
% I_yz = 0;
% 
% for i = 1:n_masses
%     x = x_CMrR_B - componentMassesAndGeom(i,5:7)';
%     m = componentMassesAndGeom(i,1);
%     dx = componentMassesAndGeom(i,2:4)';
%     I_xx = I_xx + 1/12*m*(dx(2)^2 + dx(3)^2) + m*(x(2)^2+x(3)^2);
%     I_yy = I_yy + 1/12*m*(dx(1)^2 + dx(3)^2) + m*(x(1)^2+x(3)^2);
%     I_zz = I_zz + 1/12*m*(dx(1)^2 + dx(2)^2) + m*(x(1)^2+x(2)^2);
%     I_xy = I_xy + 0 + m*x(1)*x(2);
%     I_xz = I_xz + 0 + m*x(1)*x(3);
%     I_yz = I_yz + 0 + m*x(2)*x(3);
% end
% 
% I_CM = [I_xx, I_xy, I_zz;...
%         I_xy, I_yy, I_yz;...
%         I_xz, I_yz, I_zz];
% plane.I_CM = I_CM;
% plane.mass = mass;
% plane.x_CMrR_B = x_CMrR_B;
% fprintf('Total mass of vehicle is : %8.3f [g]\n', mass)
% fprintf('Center of mass relative to reference point is located at: \n')
% fprintf('%8.3f [m]\n', x_CMrR_B')
% 
% fprintf('Inertia Tensor about center of mass is: \n')
% fprintf('[%8.3f%8.3f%8.3f]\n[%8.3f%8.3f%8.3f] [g-m^2]\n[%8.3f%8.3f%8.3f]\n', I_CM)
% % Problem 4, defining location of aero surfaces for force/moment calculation
% n_s1 = [1; 0; 0]; %dummy for now, maybe for prop later?
% n_s2 = [0; 0; -1]; 
% n_s3 = [0; 1; 0]; 
% n_s4 = [0; 0; -1]; 
% n_s5 = [0; 0; -1];
% 
% CL0_s1 = 0; %dummy for now, may be for prop later?
% CL0_s2 = 0; 
% CL0_s3 = 0; 
% CL0_s4 = 0.05; 
% CL0_s5 = 0.05;
% e_s1 = 0.9;
% e_s2 = 0.8; 
% e_s3 = 0.8; 
% e_s4 = 0.9; 
% e_s5 = 0.9;
% i_s1 = 0;
% i_s2 = 0; 
% i_s3 = 0; 
% i_s4 = 0.05; 
% i_s5 = 0.05;
% CD0_s1 = 0;
% CD0_s2 = 0.01; 
% CD0_s3 = 0.01; 
% CD0_s4 = 0.01; 
% CD0_s5 = 0.01;
% CDa_s1 = 0;
% CDa_s2 = 1; 
% CDa_s3 = 1; 
% CDa_s4 = 1; 
% CDa_s5 = 1;
% a0_s1 = 0;
% a0_s2 = 0; 
% a0_s3 = 0; 
% a0_s4 = 0.05; 
% a0_s5 = 0.05;
% CM0_s1 = 0;
% CM0_s2 = 0; 
% CM0_s3 = 0; 
% CM0_s4 = -0.05; 
% CM0_s5 = -0.05;
% CMa_s1 = 0;
% CMa_s2 = 0; 
% CMa_s3 = 0; 
% CMa_s4 = 0; 
% CMa_s5 = 0;
% s1 propeller maybe?
% c_s1 = 0;
% b_s1 = 0;
% S_s1 = 0;
% AR_s1 = 1;
% x_s1rR_B = [0;0;0];
% s2 Horizontal Stabilizer (Row 3 above)
% c_s2 = componentMassesAndGeom( 3, 2);
% b_s2 = componentMassesAndGeom( 3, 3);
% S_s2 = b_s2*c_s2;
% AR_s2 = b_s2/c_s2;
% x_s2rR_B = [componentMassesAndGeom( 3, 5) + c_s2/4;... % position of s2 lift point relative to point R in body coordinates
%             componentMassesAndGeom( 3, 6);...
%             componentMassesAndGeom( 3, 7)];
% 
% s3 Vertical Stabilizer (Row 4 above)
% c_s3 = componentMassesAndGeom( 4, 2);
% b_s3 = componentMassesAndGeom( 4, 4);   % this is the vertical stabilizer, so span is in the z dir.
% S_s3 = b_s3*c_s3;
% AR_s3 = b_s3/c_s3;
% x_s3rR_B = [componentMassesAndGeom( 4, 5) + c_s3/4;... % position of s3 lift point relative to Reference location R in body coordinates
%             componentMassesAndGeom( 4, 6);...
%             componentMassesAndGeom( 4, 7)];
% 
% s4 Right Wing (Row 1 above)
% c_s4 = componentMassesAndGeom( 1, 2);
% b_s4 = componentMassesAndGeom( 1, 3);
% S_s4 = b_s4*c_s4;
% AR_s4 = b_s4/c_s4;
% x_s4rR_B = [componentMassesAndGeom( 1, 5) + c_s4/4;... % position of s4 lift point relative to Reference location R in body coordinates
%             componentMassesAndGeom( 1, 6);...
%             componentMassesAndGeom( 1, 7)];
% 
% s5 Left Wing (Row 2 above)
% c_s5 = componentMassesAndGeom( 2, 2);
% b_s5 = componentMassesAndGeom( 2, 3);
% S_s5 = b_s5*c_s5;
% AR_s5 = b_s5/c_s5;
% x_s5rR_B = [componentMassesAndGeom( 2, 5) + c_s5/4;... % position of s5 lift point relative to Reference location R in body coordinates
%             componentMassesAndGeom( 2, 6);...
%             componentMassesAndGeom( 2, 7)];
% 
% 
% fprintf(' =======================================================================\n')
% fprintf('|         Aerodynamic Locations and Properties of Aero Surfaces         |\n')
% fprintf('|________|_________x_s#/R^B_________|___________________________________|\n')
% fprintf('|___s#___|____x___|____y___|____z___|____c___|____b___|____S___|___AR___|\n')
% fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
%     's2',x_s2rR_B, c_s2, b_s2, S_s2, AR_s2)
% fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
%     's3',x_s3rR_B, c_s3, b_s3, S_s3, AR_s3)
% fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
%     's4',x_s4rR_B, c_s4, b_s4, S_s4, AR_s4)
% fprintf('|%5s   |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |%7.3f |\n',...
%     's5',x_s5rR_B, c_s5, b_s5, S_s5, AR_s5)
% fprintf(' =======================================================================\n')
% 
% % Problem 5, compiling aerodynamic terms for use in runtime.
% position of aero surfaces WRT CM
% plane.x_srCM_B = [0;...
%             x_s2rR_B - x_CMrR_B;...
%             x_s3rR_B - x_CMrR_B;...
%             x_s4rR_B - x_CMrR_B;...
%             x_s5rR_B - x_CMrR_B;...
%             ];
% % combine all variables into single variables for ease of iteration
% plane.n   = [  n_s1,   n_s2,   n_s3,   n_s4,   n_s5];
% plane.CL0 = [CL0_s1, CL0_s2, CL0_s3, CL0_s4, CL0_s5];
% plane.e   = [  e_s1,   e_s2,   e_s3,   e_s4,   e_s5];
% plane.i_s = [  i_s1,   i_s2,   i_s3,   i_s4,   i_s5];
% plane.CD0 = [CD0_s1, CD0_s2, CD0_s3, CD0_s4, CD0_s5];
% plane.CDa = [CDa_s1, CDa_s2, CDa_s3, CDa_s4, CDa_s5];
% plane.a0  = [ a0_s1,  a0_s2,  a0_s3,  a0_s4,  a0_s5];
% plane.CM0 = [CM0_s1, CM0_s2, CM0_s3, CM0_s4, CM0_s5];
% plane.CMa = [CMa_s1, CMa_s2, CMa_s3, CMa_s4, CMa_s5];
% plane.c   = [  c_s1,   c_s2,   c_s3,   c_s4,   c_s5];
% plane.b   = [  b_s1,   b_s2,   b_s3,   b_s4,   b_s5];
% plane.S   = [  S_s1,   S_s2,   S_s3,   S_s4,   S_s5];
% plane.AR  = [ AR_s1,  AR_s2,  AR_s3,  AR_s4,  AR_s5];
% plane.n_aero = 5; % number of aero surfaces
% 
% plane.rho = 1.225;    % [kg/m^3] density of air at sea level
% 
% plane.halfRhoS = 1/2*rho*S;
% plane.halfRhoSc = halfRhoS.*c;
% 
