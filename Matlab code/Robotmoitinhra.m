clear all
clc
%L2 = 8.872; L3 = 11.5;  L4 = 4; d1 = 6.1214; d2 =3; d3 = 12.015;d4=5.26 ; theta1 = 27.0490; theta2 = -29.1991 ; theta3 = -23.4238;
%CT_JcraigRAD(alpha, a, d, theta)
syms L1 L2 L3 L4 theta1 theta2 theta3 d1 d2 d3 d4 r11 r12 r13 x r21 r22 r23 y r31 r32 r33 z
T01 = CT_JcraigRAD(0,0,d1,theta1);
T12 = CT_JcraigRAD(sym(pi/2),0,d2,theta2);
%T12 = CT_JcraigRAD(pi/2,L1,0,theta2);
T23 = CT_JcraigRAD(0,L2,d3,theta3);
T34 = CT_JcraigRAD(0,L3,d4,0);
T45 = CT_JcraigRAD(0,L4,0,0);
T05 = T01*T12*T23*T34*T45;
P = T05(:,4);
Pa = simplify(P);
theta23 = theta2 + theta3 ;
%FK = [r11 r12 r13 x; r21 r22 r23 y; r31 r32 r33 z; 0 0 0 1];
%%
Px = T05(1,4);Py= T05(2,4); Pz =T05(3,4);
%%
%Deg
%clear all
%clc
L2 = 8.872; L3 = 11.5;  L4 = 4; d1 = 6.1214; d2 =3; d3 = 12.015;d4=5.26 ; theta1 = 1.170490397160153e+02; theta2 = -28.300743682191822 ; theta3 = -21.137081153710874;
%L1 = 3;L2 = 12.15; L3 = 11.5;  L4 = 4; theta1 = 10; theta2 = 20 ; theta3 = 30 ; theta4 = 0; d1 = 6.1214; d2 =0 ; 
% = 10;L2 = 20; L3 = 25;  L4 = 5; theta1 = 10; theta2 = 29.1 ; theta3 = -28.4892 ; theta4 = -0.62; d1 = 10; d2 =0 ; 
%CT_JcraigDEG(alpha, a, d, theta)
%syms L0 L1 L2 L3 L4 theta1 theta2 theta3 theta4 d1 d2 r11 r12 r13 x r21 r22 r23 y r31 r32 r33 z
T01 = CT_JcraigDEG(0,0,d1,theta1);
%T12 = CT_JcraigRAD(sym(pi/2),L1,0,theta2);
T12 = CT_JcraigDEG(90,0,d2,theta2);
T23 = CT_JcraigDEG(0,L2,d3,theta3);
T34 = CT_JcraigDEG(0,L3,d4,0);
T45 = CT_JcraigDEG(0,L4,0,0);
T05 = T01*T12*T23*T34*T45;
P = T05(:,4);
theta23 = theta2+theta3;
%Pa = simplify(P);
%theta234 = theta2 + theta3 + theta4;
%FK = [r11 r12 r13 x; r21 r22 r23 y; r31 r32 r33 z; 0 0 0 1];
%%
Px = T05(1,4);Py= T05(2,4); Pz =T05(3,4);
%%
%DEG
%RAD
alpha =0;
alpha1 = asind(-Py/sqrt((Px^2 +Py^2)));
alpha2 = acosd(Px/sqrt((Px^2+Py^2)));
if Px>0 && Py > 0
   alpha = alpha1; 
else
   alpha = alpha2;
end
theta11 = asind((d2+d3+d4)/sqrt((Px^2 +Py^2)))-alpha;
K2 = Pz - d1-L4*sind(theta23);
K1 = Px*cosd(theta11) + Py*sind(theta11) - L4*cosd(theta23);
c3 = ((K1^2 + K2^2-L3^2-L2^2)/(2*L2*L3));
s3 = sqrt(1-(c3^2));
theta33 = -acosd(c3);
theta333 = atan2d(s3,c3);
s2 = ((K2*(L2+L3*cosd(theta33)))-(K1*L3*sind(theta33)))/((L3*sind(theta33))^2+(L2+L3*cosd(theta33))^2);
%c2 =((K1*(L2+L3*cosd(theta33)))+(K2*L3*sind(theta33)))/((L3*sind(theta33))^2+(L2+L3*cosd(theta33))^2);
c2 = (K1+L3*s2*sind(theta33))/(L2+L3*cosd(theta33));
theta22 = atan2d(s2,c2);
%theta44 = -theta2 - theta3;
%% Robot Dynamic

%%
syms L2 L3 L4  m1 m2 m3 m4 m5 theta1 theta2 theta3 d1 d2 d3 d4 alpha alpha1 alpha2 I1 I2 I3 I4 I5 fx fy fz theta1_dot theta1_ddot theta2_dot theta2_ddot theta3_dot theta3_ddot omega0 v0
assume(m4,'real'); assume(I3,'real');
assume(L2,'real');  assume(m1,'real');  assume(m5,'real'); assume(I4,'real');
assume(L3,'real');  assume(m2,'real');  assume(I1,'real'); assume(I5,'real');
assume(L4,'real');  assume(m3,'real');  assume(I2,'real');
assume(theta1,'real'); assume(theta1_dot,'real'); assume(theta1_ddot,'real');
assume(theta2,'real'); assume(theta2_dot,'real'); assume(theta2_ddot,'real');
assume(theta3,'real'); assume(theta3_dot,'real'); assume(theta3_ddot,'real');
assume(d1,'real');
assume(d2,'real');
assume(d3,'real');
assume(d4,'real');
assume(alpha,'real');
assume(alpha1,'real');
assume(alpha2,'real');
assume(fx,'real');
assume(fy,'real');
assume(fz,'real');
assume(omega0,'real');
R01 = T01(1:3,1:3); R011 = R01'; P01 = T01(1:3,4); 
R12 = T12(1:3,1:3); R122 = R12'; P12 = T12(1:3,4);
R23 = T23(1:3,1:3); R233 = R23'; P23 = T23(1:3,4); 
R34 = T34(1:3,1:3); R344 = R34'; P34 = T34(1:3,4);
R45 = T45(1:3,1:3); R455 = R45'; P45 = T45(1:3,4);
R02 = R01*R12;   T02 = T01*T12;
R03 = R02*R23;   T03 = T02*T23;
R04 = R03*R34;   T04 = T03*T34;
R05 = R04*R45;   T05 = T04*T45;
%%
%% LINK VELOCITy
omega0 = [0;0;0]; v0 = [0;0;0];
P0c = [0 0 d1]'; P1c = [0 -d2 0]'; P2c = [L2/2 0 d3]'; P3c = [L3/2 0 d4]'; P4c = [L4/2 0 0]';
theta_dot = [theta1_dot;theta2_dot;theta3_dot;0;0];
R = {R011,R122, R233, R344, R455};
Rc = {R01,R02,R03,R04,R05};
P = {P01,P12,P23,P34,P45};
Pc = {P0c,P1c,P2c,P3c,P4c};
T = {T01,T02,T03,T04,T05};
omega = sym(zeros(3, 5));
v = sym(zeros(3, 5));
vc = sym(zeros(3,5));
Pcc = sym(zeros(4,5));
for i = 0:4
    switch i
        case 0
            thetadoti = [0;0;theta_dot(1)];
        case 1
            thetadoti = [0;0;theta_dot(2)];
        case 2
            thetadoti = [0;0;theta_dot(3)];
        case 3
            thetadoti = [0;0;theta_dot(4)];
        case 4
            thetadoti = [0;0;theta_dot(5)];
    end
        if i == 0
            omega(:, i+1) = R{i+1} * omega0 +thetadoti;
            v(:, i+1) = R{i+1} * v0 + SKMN(omega0) * P{i+1};
            vc(:,i+1) = Rc{i+1} * v(:,i+1) + SKMN(omega(:,i+1)) * Pc{i+1};
        else
            omega(:, i+1) = R{i+1} * omega(:, i)+thetadoti;
            v(:, i+1) = R{i+1} * v(:, i) + SKMN(omega(:, i)) * P{i+1};
            vc(:,i+1) = Rc{i+1}*v(:,i+1)+SKMN(omega(:,i+1))*Pc{:,i+1};
        end
        omega(:, i+1) = simplify(omega(:, i+1));
        v(:, i+1) = simplify(v(:, i+1));
        vc(:,i+1) = simplify(vc(:,i+1));
end 
 Pcc(:,1)= T{1}*[P0c;1];
 Pcc(:,2)= T{2}*[P1c;1];
 Pcc(:,3)= T{3}*[P2c;1];
 Pcc(:,4)= T{4}*[P3c;1];
 Pcc(:,5)= T{5}*[P4c;1];
 Pcc = simplify(Pcc);
%% Larange
syms g t;

K = sym(zeros(1, 5));
U = sym(zeros(1, 5));
m = [m1, m2, m3, m4, m5];
I = [I1, I2, I3, I4, I5];
g = [0, g, 0];

for i = 1:5
    % Tính n?ng l??ng cinetic c?a m?i ph?n t?
    K(i) = 1/2 * m(i) * vc(:,i)' * vc(:,i) + I(i) * omega(:,i)' * omega(:,i);
    
    % Tính n?ng l??ng ti?m n?ng c?a m?i ph?n t?
    U(i) = -m(i) * g * Pcc(1:3,i);
end

% T?ng n?ng l??ng cinetic và n?ng l??ng ti?m n?ng c?a toàn b? h? th?ng
Ktotal = sum(K);
Utotal = sum(U);

% Tính hàm Lagrangian
L = simplify(Ktotal - Utotal);

%%
% Define symbolic variables for theta_dot and theta_ddot
theta_dot = [theta1_dot; theta2_dot; theta3_dot];
theta_ddot = [theta1_ddot; theta2_ddot; theta3_ddot];
theta = [theta1; theta2; theta3];

% Calculate the derivatives of the Lagrangian with respect to theta_dot and theta
dLdtheta_dot = jacobian(L, theta_dot).';
dLdtheta = jacobian(L, theta).';

% Calculate the time derivatives of dLdtheta_dot
dLdtheta_dot_dt = jacobian(dLdtheta_dot, theta) * theta_dot + jacobian(dLdtheta_dot, theta_dot) * theta_ddot;

% Calculate the equations of motion
tau = simplify(dLdtheta_dot_dt - dLdtheta);

% Extract individual torques
te11 = tau(1);
te22 = tau(2);
te33 = tau(3);
%%
% Extract the inertia matrix M
M = jacobian(tau, theta_ddot);
M = simplify(M);

% Extract the Coriolis and centrifugal matrix V
C = sym(zeros(3, 3));
for k = 1:3
    for j = 1:3
        for i = 1:3
            C(k, j) = C(k, j) + 0.5 * (diff(M(k, j), theta(i)) + diff(M(k, i), theta(j)) - diff(M(i, j), theta(k))) * theta_dot(i);
        end
    end
end
V = simplify(C * theta_dot);

% Extract the gravity vector G
G = jacobian(Utotal, theta).';
G = simplify(G);

% Display the results
disp('Inertia Matrix (M):');
disp(M);

disp('Coriolis and Centrifugal Matrix (V):');
disp(V);

disp('Gravity Vector (G):');
disp(G);

%%