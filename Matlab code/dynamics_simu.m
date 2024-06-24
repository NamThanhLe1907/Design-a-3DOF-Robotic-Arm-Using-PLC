function [M, V, G] = dynamics_simu(theta, theta_dot, theta_ddot, L2, L3, L4, d1, d2, d3, d4, mm, II, g)
    % Kinematic and dynamic parameters
    theta1 = theta(1); theta2 = theta(2); theta3 = theta(3);
    theta1_dot = theta_dot(1); theta2_dot = theta_dot(2); theta3_dot = theta_dot(3);
    theta1_ddot = theta_ddot(1); theta2_ddot = theta_ddot(2); theta3_ddot = theta_ddot(3);
    m1 = mm(1); m2 =mm(2); m3 = mm(3); m4 = mm(4); m5 = mm(5);
    I1 = II(1); I2 = II(2); I3 = II(3); I4 = II(4); I5 = II(5);
    % Define necessary matrices and vectors
    K = zeros(1, 5);
    U = zeros(1, 5);
    omega = zeros(3, 5);
    v = zeros(3, 5);
    vc = zeros(3,5);
    omega0 = [0;0;0]; v0 = [0;0;0];
    P0c = [0 0 d1]'; P1c = [0 -d2 0]'; P2c = [L2/2 0 d3]'; P3c = [L3/2 0 d4]'; P4c = [L4/2 0 0]';
    theta_dot = [theta1_dot;theta2_dot;theta3_dot;0;0];
    R = {R011,R122, R233, R344, R455};
    Rc = {R01,R02,R03,R04,R05};
    P = {P01,P12,P23,P34,P45};
    Pc = {P0c,P1c,P2c,P3c,P4c};
    T = {T01,T02,T03,T04,T05};
    Pcc = zeros(4,5);
    m = [m1, m2, m3, m4, m5];
    I = [I1, I2, I3, I4, I5];
    g = [0,9.8,0];
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
    end 
    Pcc(:,1)= T{1}*[P0c;1];
    Pcc(:,2)= T{2}*[P1c;1];
    Pcc(:,3)= T{3}*[P2c;1];
    Pcc(:,4)= T{4}*[P3c;1];
    Pcc(:,5)= T{5}*[P4c;1];
    % Compute Lagrangian and equations of motion as previously shown
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
    L = Ktotal - Utotal;
    % Calculate the derivatives of the Lagrangian
    % Calculate the derivatives of the Lagrangian with respect to theta_dot and theta
    dLdtheta_dot = jacobian(L, theta_dot).';
    dLdtheta = jacobian(L, theta).';

    % Calculate the time derivatives of dLdtheta_dot
    dLdtheta_dot_dt = jacobian(dLdtheta_dot, theta) * theta_dot + jacobian(dLdtheta_dot, theta_dot) * theta_ddot;
    tau = dLdtheta_dot_dt - dLdtheta;
    % Extract individual torques and matrices
    M = jacobian(tau, theta_ddot);
    C = zeros(3, 3);
    for k = 1:3
        for j = 1:3
            for i = 1:3
                C(k, j) = C(k, j) + 0.5 * (diff(M(k, j), theta(i)) + diff(M(k, i), theta(j)) - diff(M(i, j), theta(k))) * theta_dot(i);
            end
        end
    end
    V = C * theta_dot;
    G = jacobian(Utotal, theta).';
end
