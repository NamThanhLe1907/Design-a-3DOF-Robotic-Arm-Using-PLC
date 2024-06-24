%% Robot dynamic
clear all
clc
%L1 = 10;L2 = 20; L3 = 25;  L4 = 5; theta1 = 0.1; theta2 = 0.2 ; theta3 = -0.3 ; 0 = 0.1; d1 = 10; d2 =0 ; 
%CT_JcraigRAD(alpha, a, d, theta)
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
%clear all
%clc
%L1 = 10;L2 = 20; L3 = 25;  L4 = 5; theta1 = 0.1; theta2 = 0.756 ; theta3 = -0.2 ; 0 = 0.4433; d1 = 10; d2 =0 ; 
%L1 = 10;L2 = 20; L3 = 25;  L4 = 5; theta1 = 0.1; theta2 = 0.1914 ; theta3 = -0.2846 ; 0 = 0.0932; d1 = 10; d2 =0 ; 
%CT_JcraigRAD(alpha, a, d, theta)
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
Px = T05(1,4);Py= T05(2,4); Pz =T05(3,4);
%%
R01 = T01(1:3,1:3); R011 = R01'; P01 = T01(1:3,4); 
R12 = T12(1:3,1:3); R122 = R12'; P12 = T12(1:3,4);
R23 = T23(1:3,1:3); R233 = R23'; P23 = T23(1:3,4); 
R34 = T34(1:3,1:3); R344 = R34'; P34 = T34(1:3,4);
R45 = T45(1:3,1:3); R455 = R45'; P45 = T45(1:3,4);
R02 = R01*R12;   T02 = T01*T12;
R03 = R02*R23;   T03 = T02*T23;
R04 = R03*R34;   T04 = T03*T34;
R05 = R04*R45;   T05 = T04*T45;
%% Robot Dynamic Static Force and Torque
f6 = [fx;fy;fz];
n6 = [0;0;0];
P5=SKM(L4,0,0);
P4=SKM(L3,0,d4);
P3=SKM(L2,0,d3); 
P2=SKM(0,-d2,0);
P1=SKM(0,0,d1);
matrixZ = [0;0;1]; 
R = {R01,R12, R23, R34, R45};
P = {P1,P2,P3,P4,P5};
f = sym(zeros(3, 5));
n = sym(zeros(3, 5));
tn = sym(zeros(5, 1));
tf = sym(zeros(5, 1));
for i = 5:-1:1
    if i == 5
        fi = fa(R{i},f6);
        ni = na(R{i},n6,P{i},fi);
    else
        fi = fa(R{i}, f(:, i+1));
        ni = na(R{i}, n(:, i+1), P{i}, fi);
    end
    f(:, i) = fi;
    n(:, i) = ni;
    % Assuming f() and n() are functions that take Ri, fi, and ni as inputs
    tn(i) = n(:, i)' * matrixZ;
    tf(i) = f(:, i)' * matrixZ;
end 
f = simplify(f);
n = simplify(n);
tn = simplify(tn);
tf = simplify(tf);
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
%% JABCOBIAN
R05 = R01*R12*R23*R34*R45; R05 = simplify(R05); f05 = R05*f(:,4);f05 = simplify(f05);
Jv5 = [L1*sin(theta3+0), L2*sin(0), 0 ,0;
       L1*cos(theta3+0), L3+L4+L2*cos(0), L3+L4, L4;
       -(L3*cos(theta2+theta3)+L2*cos(theta2)+L4*cos(theta2+ theta3)),0,0,0];
Jv0 = simplify(R05*Jv5);
tor = Jv0'*f05;
tor = simplify(tor);
tor1 = Jv5'*f(:,4);
tor1 = simplify(tor1);
%% NEWTON-EULER
omega_dot0 = [0;0;0];
omega_dot1 = R011*(omega_dot0+cross(omega0,[0;0;theta1_dot])+[0;0;theta1_ddot]);        omega_dot1 = simplify(omega_dot1);
omega_dot2 = R122*(omega_dot1+cross(omega(:,2),[0;0;theta1_dot])+[0;0;theta1_ddot]);    omega_dot2 = simplify(omega_dot2);
omega_dot3 = R233*(omega_dot2+cross(omega(:,3),[0;0;theta2_dot])+[0;0;theta2_ddot]);    omega_dot3 = simplify(omega_dot3);
omega_dot4 = R344*(omega_dot3+cross(omega(:,4),[0;0;theta3_dot])+[0;0;theta3_ddot]);    omega_dot4 = simplify(omega_dot4);
omega_dot5 = R455*(omega_dot4+cross(omega(:,5),[0;0;0])+[0;0;0]);    omega_dot5 = simplify(omega_dot5);

%% Linear acceleration
a0 = [0;g;0];  
P0c = [0;0;0]; P1c = [L1/2 0 0]'; P2c = [L2/2 0 0]'; P3c = [L3/2 0 0]'; P4c = [L4/2 0 0]'; 
a1 = R011*(a0+cross(omega_dot0,P01)+cross(omega_dot0,cross(omega_dot0,P01))); a1 = simplify(a1);
a2 = R122*(a1+cross(omega_dot1,P12)+cross(omega_dot1,cross(omega_dot1,P12))); a2 = simplify(a2);
a3 = R233*(a2+cross(omega_dot2,P23)+cross(omega_dot2,cross(omega_dot2,P23))); a3 = simplify(a3);
a4 = R344*(a3+cross(omega_dot3,P34)+cross(omega_dot3,cross(omega_dot3,P34))); a4 = simplify(a4);
a5 = R455*(a4+cross(omega_dot4,P45)+cross(omega_dot4,cross(omega_dot4,P45))); a5 = simplify(a5);   

%% Linear acceleration at center-point
ac1 = a1 + cross(omega_dot1,P0c) +cross(omega_dot1,cross(omega_dot1,P0c)); ac1 = simplify(ac1);
ac2 = a2 + cross(omega_dot2,P1c) +cross(omega_dot2,cross(omega_dot2,P1c)); ac2 = simplify(ac2);
ac3 = a3 + cross(omega_dot3,P2c) +cross(omega_dot3,cross(omega_dot3,P2c)); ac3 = simplify(ac3);
ac4 = a4 + cross(omega_dot4,P3c) +cross(omega_dot4,cross(omega_dot4,P3c)); ac4 = simplify(ac4);
ac5 = a5 + cross(omega_dot5,P4c) +cross(omega_dot5,cross(omega_dot5,P4c)); ac5 = simplify(ac5);
%% Force-balance and Joint Force
F1 = m1*ac1; f1 = f(:,4)+F1; f1 = simplify(f1);
F2 = m2*ac2; f2 = f(:,3)+F2; f2 = simplify(f2);
F3 = m3*ac3; f3 = f(:,2)+F3; f3 = simplify(f3);
F4 = m4*ac4; f4 = f(:,1)+F4; f4 = simplify(f4);
F5 = m5*ac5; f5 = f5+F4; f5 = simplify(f5);
%% Torque-Balance 
Ic1 = 0; Ic2 = 0; Ic3 = 0; Ic4 = 0; Ic5 = 0; %Point-mass assumption
N1 = Ic1*omega_dot1+cross(omega_dot1,Ic1*omega_dot1);
N2 = Ic2*omega_dot2+cross(omega_dot2,Ic2*omega_dot2);
N3 = Ic3*omega_dot3+cross(omega_dot3,Ic3*omega_dot3);
N4 = Ic4*omega_dot4+cross(omega_dot4,Ic4*omega_dot4);
N5 = Ic5*omega_dot5+cross(omega_dot5,Ic5*omega_dot5);
%% Torque equation
nn5= N5 + R45*n5+cross(P4c,F5)+cross(P45,R45*f5); nn5 = simplify(nn5);
nn4= N4 + R34*nn5+cross(P3c,F4)+cross(P34,R34*f4); nn4 = simplify(nn4);
nn3= N3 + R23*nn4+cross(P2c,F3)+cross(P23,R23*f3); nn3 = simplify(nn3);
nn2= N2 + R12*nn3+cross(P1c,F2)+cross(P12,R12*f2); nn2 = simplify(nn2);
nn1= N1 + R01*nn2+cross(P0c,F1)+cross(P01,R01*f1); nn1 = simplify(nn1);
%% New
% Khai báo d? li?u
syms omega_dot a ac F N g 
assume(omega_dot,'real');
assume(a,'real');
assume(g,'real');
assume(ac,'real');
assume(F,'real');
assume(N,'real');
m = [m1, m2, m3, m4, m5];
omega_dot0 = [0;0;0];
omega_dot = sym(zeros(3, 6)); 
a0 = [0; g; 0];  
P0c = [0 0 d1]'; P1c = [0 -d2 0]'; P2c = [L2/2 0 d3]'; P3c = [L3/2 0 d4]'; P4c = [L4/2 0 0]';
P = {P01,P12,P23,P34,P45};
a = sym(zeros(3, 6));
ac = sym(zeros(3, 6));
F = sym(zeros(3, 6));
N = sym(zeros(3, 6));
Pc = {P0c, P1c, P2c, P3c, P4c};
theta_ddot = [theta1_ddot; theta2_ddot; theta3_ddot; 0; 0];
ff = {f(:,5),f(:,4),f(:,3),f(:,2),f(:,1)};
nn = sym(zeros(3, 6)); % Kh?i t?o ma tr?n nn
for i = 1:5
    thetadoti = [0; 0; theta_dot(i)];
    thetaddoti = [0; 0; theta_ddot(i)];
    
    if i == 1
        omega_dot(:, i) = R011*(omega_dot0 + cross(omega0,[0;0;theta1_dot]) + [0;0;theta1_ddot]);
        a(:, i) = R011*(a0 + cross(omega_dot0,P01) + cross(omega_dot0,cross(omega_dot0,P01)));
    else
        omega_dot(:, i) = R{i}*(omega_dot(:, i-1) + cross(omega(:,i),thetadoti) + thetaddoti);
        a(:, i) = R{i}*(a(:, i-1) + cross(omega_dot(:, i-1),P{i}) + cross(omega_dot(:, i-1),cross(omega_dot(:, i-1),P{i})));
    end
    omega_dot(:, i) = simplify(omega_dot(:, i));
    a(:, i) = simplify(a(:, i));
    
    % Linear acceleration at center-point
    ac(:, i) = a(:, i) + cross(omega_dot(:, i), Pc{i}) + cross(omega_dot(:, i), cross(omega_dot(:, i), Pc{i}));
    ac(:, i) = simplify(ac(:, i));
    
    % Force-balance and Joint Force
    F(:, i) = m(i) * ac(:, i);
    ff{:, i} = ff{:, 6 - i} + F(:, i);
    ff{:, i} = simplify(ff{:, i});
    
    % Torque-Balance 
    Ic = 0; % Point-mass assumption
    N(:, i) = Ic * omega_dot(:, i) + cross(omega_dot(:, i), Ic * omega_dot(:, i));
    N(:, i) = simplify(N(:, i));
    %Torque
end
%%
    for j = 5:-1:1
    if j == 5
    nn(:, j) = N(:, 5) + R{5}*n6 + cross(Pc{5}, F(:, 5)) + cross(P{5}, R{5}*f6); 
    nn(:, j) = simplify(nn(:, j));
    end
    nn(:, j) = N(:, j) + R{j}*nn(:, j+1) + cross(Pc{j}, F(:, j)) + cross(P{j}, R{j}*ff{:, j});
    nn(:, j) = simplify(nn(:, j));
    end
tau1 = nn(:, 5);
tau2 = nn(:, 4);
tau3 = nn(:, 3);
tau4 = nn(:, 2);

%M = equationsToMatrix([tau1; tau2; tau3; tau4], [theta1_ddot; theta2_ddot; theta3_ddot; 0]);
%G = equationsToMatrix([tau1; tau2; tau3; tau4], [g]);
%%
tau1 = nn(:,5);
tau2 = nn(:,4);
tau3 = nn(:,3);
tau4 = nn(:,2);
tau = [tau1, tau2, tau3, tau4]; % Corrected definition

% Define the vector of second-order accelerations
syms theta1_ddot theta2_ddot theta3_ddot 0

theta_ddot = [theta1_ddot; theta2_ddot; theta3_ddot; 0];

% Compute the mass matrix M
M = sym(zeros(12)); % Initialize the mass matrix

for i = 1:4
    for j = 1:4
        % Compute the derivative of the torque with respect to the j-th acceleration
        dtau_dtheta_ddot = simplify(jacobian(tau(:,i), theta_ddot(j)));
            switch i
                case 1
                    M(1:3,j) = dtau_dtheta_ddot;
                case 2
                    M(4:6,j) = dtau_dtheta_ddot;
                case 3
                    M(7:9,j) = dtau_dtheta_ddot;
                case 4
                    M(10:12,j) = dtau_dtheta_ddot;
            end

    end
end

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
% Tính các ??o hàm c?n thi?t
dLddt1 = diff(L, theta1_dot);
dLddt2 = diff(L, theta2_dot);
dLddt3 = diff(L, theta3_dot);

dLdt1 = diff(L, theta1);
dLdt2 = diff(L, theta2);
dLdt3 = diff(L, theta3);
%%
dLddt1_dt = diff(dLddt1, theta1) * theta1_dot + diff(dLddt1, theta2) * theta2_dot + diff(dLddt1, theta3) * theta3_dot ...
    + diff(dLddt1, theta1_dot) * theta1_ddot + diff(dLddt1, theta2_dot) * theta2_ddot + diff(dLddt1, theta3_dot) * theta3_ddot;

dLddt2_dt = diff(dLddt2, theta1) * theta1_dot + diff(dLddt2, theta2) * theta2_dot + diff(dLddt2, theta3) * theta3_dot ...
    + diff(dLddt2, theta1_dot) * theta1_ddot + diff(dLddt2, theta2_dot) * theta2_ddot + diff(dLddt2, theta3_dot) * theta3_ddot;

dLddt3_dt = diff(dLddt3, theta1) * theta1_dot + diff(dLddt3, theta2) * theta2_dot + diff(dLddt3, theta3) * theta3_dot ...
    + diff(dLddt3, theta1_dot) * theta1_ddot + diff(dLddt3, theta2_dot) * theta2_ddot + diff(dLddt3, theta3_dot) * theta3_ddot;

%%
syms te1 te2 te3 te4 t
% Tính l?c theo ph??ng pháp Lagrange-Euler
te1 = simplify(dLddt1_dt - dLdt1);
te2 = simplify(dLddt2_dt - dLdt2);
te3 = simplify(dLddt3_dt - dLdt3);


% Rút g?n các bi?u th?c
te1 = simplify(te1);
te2 = simplify(te2);
te3 = simplify(te3);
%% Larange
t1 = simplify(diff(dLddt1,theta1)*theta1_dot+diff(dLddt1,theta2)*theta2_dot+diff(dLddt1,theta3)*theta3_dot ...
    +diff(dLddt1,theta1_dot)*theta1_ddot+diff(dLddt1,theta2_dot)*theta2_ddot+diff(dLddt1,theta3_dot)*theta3_ddot) ;
t2 = simplify(diff(dLddt2,theta1)*theta1_dot+diff(dLddt2,theta2)*theta2_dot+diff(dLddt2,theta3)*theta3_dot ...
    +diff(dLddt2,theta1_dot)*theta1_ddot+diff(dLddt2,theta2_dot)*theta2_ddot+diff(dLddt2,theta3_dot)*theta3_ddot) ;
t3 = simplify(diff(dLddt3,theta1)*theta1_dot+diff(dLddt3,theta2)*theta2_dot+diff(dLddt3,theta3)*theta3_dot ...
    +diff(dLddt3,theta1_dot)*theta1_ddot+diff(dLddt3,theta2_dot)*theta2_ddot+diff(dLddt3,theta3_dot)*theta3_ddot) ;
%% M matrix
M = simplify([diff(t1,theta1_ddot),diff(t1,theta2_ddot),diff(t1,theta3_ddot);
    diff(t2,theta1_ddot),diff(t2,theta2_ddot),diff(t2,theta3_ddot);
    diff(t3,theta1_ddot),diff(t3,theta2_ddot),diff(t3,theta3_ddot);]);
%% G Matrix
G = simplify([diff(Utotal,theta1);...
    diff(Utotal,theta2);...
    diff(Utotal,theta3);...
    diff(Utotal,0)]);
%% Coriolis Matrix
dM = simplify([diff(M(1,1),theta1)*theta1_dot+diff(M(1,1),theta2)*theta2_dot+diff(M(1,1),theta3)*theta3_dot+diff(M(1,1),0)*0+diff(M(1,2),theta1)*theta1_dot+diff(M(1,2),theta2)*theta2_dot+diff(M(1,2),theta3)*theta3_dot+diff(M(1,2),0)*0+diff(M(1,3),theta1)*theta1_dot+diff(M(1,3),theta2)*theta2_dot+diff(M(1,3),theta3)*theta3_dot+diff(M(1,3),0)*0++diff(M(1,4),theta1)*theta1_dot+diff(M(1,4),theta2)*theta2_dot+diff(M(1,4),theta3)*theta3_dot+diff(M(1,4),0)*0;...
    diff(M(2,1),theta1)*theta1_dot+diff(M(2,1),theta2)*theta2_dot+diff(M(2,1),theta3)*theta3_dot+diff(M(2,1),0)*0+diff(M(2,2),theta1)*theta1_dot+diff(M(2,2),theta2)*theta2_dot+diff(M(2,2),theta3)*theta3_dot+diff(M(2,2),0)*0+diff(M(2,3),theta1)*theta1_dot+diff(M(2,3),theta2)*theta2_dot+diff(M(2,3),theta3)*theta3_dot+diff(M(2,3),0)*0++diff(M(2,4),theta1)*theta1_dot+diff(M(2,4),theta2)*theta2_dot+diff(M(2,4),theta3)*theta3_dot+diff(M(2,4),0)*0;...
diff(M(3,1),theta1)*theta1_dot+diff(M(3,1),theta2)*theta2_dot+diff(M(3,1),theta3)*theta3_dot+diff(M(3,1),0)*0+diff(M(3,2),theta1)*theta1_dot+diff(M(3,2),theta2)*theta2_dot+diff(M(3,2),theta3)*theta3_dot+diff(M(3,2),0)*0+diff(M(3,3),theta1)*theta1_dot+diff(M(3,3),theta2)*theta2_dot+diff(M(3,3),theta3)*theta3_dot+diff(M(3,3),0)*0++diff(M(3,4),theta1)*theta1_dot+diff(M(3,4),theta2)*theta2_dot+diff(M(3,4),theta3)*theta3_dot+diff(M(3,4),0)*0;...
diff(M(4,1),theta1)*theta1_dot+diff(M(4,1),theta2)*theta2_dot+diff(M(4,1),theta3)*theta3_dot+diff(M(4,1),0)*0+diff(M(4,2),theta1)*theta1_dot+diff(M(4,2),theta2)*theta2_dot+diff(M(4,2),theta3)*theta3_dot+diff(M(4,2),0)*0+diff(M(4,3),theta1)*theta1_dot+diff(M(4,3),theta2)*theta2_dot+diff(M(4,3),theta3)*theta3_dot+diff(M(4,3),0)*0++diff(M(4,4),theta1)*theta1_dot+diff(M(4,4),theta2)*theta2_dot+diff(M(4,4),theta3)*theta3_dot+diff(M(4,4),0)*0]);

dMdt1 = Comp_dM(M,theta1); 
dMdt2 = Comp_dM(M,theta2);
dMdt3 = Comp_dM(M,theta3);
theta_dot = [theta1_dot,theta2_dot,theta3_dot,0]';
theta_ddot = [theta1_ddot,theta2_ddot,theta3_ddot,0]';
t = [t1,t2,t3,t4]';
C = simplify(dM-1/2*([theta_dot'*dMdt1;theta_dot'*dMdt2;theta_dot'*dMdt3;theta_dot'*dMdt4;]));
%%
tee = simplify(t-M*theta_ddot-C*theta_dot-G);
teee = simplify(M*theta_ddot+C*theta_dot+G);
%%
dkddt1 = simplify(diff(Ktotal,theta1_dot));
dkddt2 = simplify(diff(Ktotal,theta2_dot));
dkddt3 = simplify(diff(Ktotal,theta3_dot));
dkddt4 = simplify(diff(Ktotal,0));
DK = [dkddt1;dkddt2;dkddt3;dkddt4];
dkdt1 = simplify(diff(Ktotal,theta1));
dkdt2 = simplify(diff(Ktotal,theta2));
dkdt3 = simplify(diff(Ktotal,theta3));
dkdt4 = simplify(diff(Ktotal,0));
DKT = [dkdt1;dkdt2;dkdt3;dkdt4];

dudt1 = simplify(diff(Utotal,theta1));
dudt2 = simplify(diff(Utotal,theta2));
dudt3 = simplify(diff(Utotal,theta3));
dudt4 = simplify(diff(Utotal,0));
DU = [dudt1;dudt2;dudt3;dudt4];
