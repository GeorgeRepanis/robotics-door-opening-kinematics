% Εργασία Ρομποτική 2025 - Προσομοίωση πόρτας με Unit Quaternion καταγραφή (συμβατή με όλες τις εκδόσεις)

clc;
clear;
close all;

% Χρονικές παράμετροι
T = 5;                  % συνολικός χρόνος (sec)
dt = 0.05;              % βήμα
t = 0:dt:T;             % χρονική διακριτοποίηση

% Καθορισμός ορίων σταδίων
T1 = T/3;
T2 = 2*T/3;

% Quintic profiles για τα στάδια
s1 = 10*(t(t<=T1)/T1).^3 - 15*(t(t<=T1)/T1).^4 + 6*(t(t<=T1)/T1).^5;  
s2 = 10*((t(t>T1 & t<=T2)-T1)/(T2-T1)).^3 - 15*((t(t>T1 & t<=T2)-T1)/(T2-T1)).^4 + 6*((t(t>T1 & t<=T2)-T1)/(T2-T1)).^5; 
s3 = 10*((t(t>T2)-T2)/(T-T2)).^3 - 15*((t(t>T2)-T2)/(T-T2)).^4 + 6*((t(t>T2)-T2)/(T-T2)).^5; 

% Διαστάσεις
l = 1;                  
h_porta = 2;            
h = 0.7;                
lo = 0.1;               

% Θέση πλαισίου H ως προς D
p_dh = [l - lo; 0; h];   

% Αποθήκευση τροχιάς
traj_pos = zeros(3, length(t));    % θέσεις του H
traj_quat = zeros(4, length(t));   % quaternions του H

% Ρύθμιση γραφικών
figure(1);
axis equal;
grid on;
view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-1.5 2]);
ylim([-1.5 2]);
zlim([0 2.5]);
title('Προσομοίωση ανοίγματος πόρτας και σωστής κίνησης πλαισίου {H}');
hold on;

% Ορισμός σταθερών
T_d = transl(0,0,0); 
[X, Z] = meshgrid([0 l], [0 h_porta]); 

% Μετρητές για στάδια
i1 = 1; i2 = 1; i3 = 1;

% Προσομοίωση
for k = 1:length(t)
    
    theta = 0; phi = 0; corr = 0;
    
    if t(k) <= T1
        phi = -pi/4 * s1(i1); 
        i1 = i1 + 1;
    elseif t(k) <= T2
        phi = -pi/4;           
        theta = -pi/6 * s2(i2); 
        i2 = i2 + 1;
    else
        theta = -pi/6; 
        phi = -pi/4;   
        corr = pi/6 * s3(i3); 
        i3 = i3 + 1;
    end
    
    % --- Υπολογισμός συνολικού μετασχηματισμού ---
    Rz = trotz(theta);
    Rz3 = Rz(1:3,1:3);
    
    Rx = trotx(phi);
    Rx3 = Rx(1:3,1:3);
    
    Rcorr = trotz(corr);
    Rcorr3 = Rcorr(1:3,1:3);
    
    R_dh = Rz3 * Rx3 * Rcorr3;
    
    p = Rz3 * p_dh;
    
    T_dh = [R_dh, p; 0 0 0 1];
    
    % ---- Καταγραφή τροχιάς ----
    traj_pos(:,k) = p; % θέση
    q = rotm2quat_custom(R_dh); % υπολογισμός quaternion
    traj_quat(:,k) = q(:);
    
    % ---- Σχεδίαση πόρτας ----
    porta_vertices = [0 0 0;
                      l 0 0;
                      l 0 h_porta;
                      0 0 h_porta]';
    porta_transformed = Rz3 * porta_vertices;
    plot3([porta_transformed(1,:) porta_transformed(1,1)], ...
          [porta_transformed(2,:) porta_transformed(2,1)], ...
          [porta_transformed(3,:) porta_transformed(3,1)], ...
          'k-', 'LineWidth', 2);
    
    % ---- Σχεδίαση πλαισίου {H} ----
    trplot(T_dh, 'frame', 'H', 'color', 'r', 'length', 0.3);
    
    % ---- Πλαίσια {0} και {D} ----
    trplot(eye(4), 'frame', '0', 'color', 'k', 'length', 0.5);
    trplot(T_d, 'frame', 'D', 'color', 'b', 'length', 0.4);
    
    pause(0.01);
    
    if k < length(t)
        cla;
    end
end

% ---- Σχεδίαση τροχιάς θέσης ----
figure(2);
plot3(traj_pos(1,:), traj_pos(2,:), traj_pos(3,:), 'r', 'LineWidth', 2);
grid on;
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Τροχιά Θέσης του Πλαισίου {H}');

% ---- Σχεδίαση quaternion συνιστωσών ----
figure(3);
plot(t, traj_quat(1,:), 'r', t, traj_quat(2,:), 'g', t, traj_quat(3,:), 'b', t, traj_quat(4,:), 'k', 'LineWidth', 2);
grid on;
xlabel('Χρόνος (sec)');
ylabel('Quaternion Συνιστώσα');
legend('q_0', 'q_1', 'q_2', 'q_3');
title('Εξέλιξη Τεταρτογωνίων (Unit Quaternions) του Πλαισίου {H}');

% ---- Υποβοηθητική συνάρτηση: Υπολογισμός quaternion από rotation matrix ----
function q = rotm2quat_custom(R)
% Υπολογισμός Unit Quaternion από 3x3 rotation matrix R

tr = trace(R);
if tr > 0
    S = sqrt(tr+1.0) * 2;
    qw = 0.25 * S;
    qx = (R(3,2) - R(2,3)) / S;
    qy = (R(1,3) - R(3,1)) / S;
    qz = (R(2,1) - R(1,2)) / S;
else
    if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2;
        qw = (R(3,2) - R(2,3)) / S;
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S;
        qz = (R(1,3) + R(3,1)) / S;
    elseif (R(2,2) > R(3,3))
        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2;
        qw = (R(1,3) - R(3,1)) / S;
        qx = (R(1,2) + R(2,1)) / S;
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2)) / S;
    else
        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2;
        qw = (R(2,1) - R(1,2)) / S;
        qx = (R(1,3) + R(3,1)) / S;
        qy = (R(2,3) + R(3,2)) / S;
        qz = 0.25 * S;
    end
end
q = [qw qx qy qz];
end
