function [T01,T12,T23,T34,T45,T05] = forward_kinematics_simu(theta1, theta2, theta3, d1, d2, d3, L2, L3, L4)
    % Define the transformation matrices using Craig's convention
    T01 = CT_JcraigDEG(0, 0, d1, theta1);
    T12 = CT_JcraigDEG(90, 0, d2, theta2);
    T23 = CT_JcraigDEG(0, L2, d3, theta3);
    T34 = CT_JcraigDEG(0, L3, d4, 0);
    T45 = CT_JcraigDEG(0, L4, 0, 0);
    R01 = T01(1:3,1:3); R011 = R01'; P01 = T01(1:3,4); 
    R12 = T12(1:3,1:3); R122 = R12'; P12 = T12(1:3,4);
    R23 = T23(1:3,1:3); R233 = R23'; P23 = T23(1:3,4); 
    R34 = T34(1:3,1:3); R344 = R34'; P34 = T34(1:3,4);
    R45= T45(1:3,1:3); R455 = R45'; P45 = T45(1:3,4);
    R02 = R01*R12;   T02 = T01*T12;
    R03 = R02*R23;   T03 = T02*T23;
    R04 = R03*R34;   T04 = T03*T34;
    R05 = R04*R45;   T05 = T04*T45;
    % Compute the overall transformation matrix
    T05 = T01 * T12 * T23 * T34 * T45;
end
