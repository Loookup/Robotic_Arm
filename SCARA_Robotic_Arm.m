%% Initialize

close all;
clc;
clear;


%% Symbolic parameter

% theta

syms th1 th2;

% theta_prime

syms th_p1 th_p2;

% Theta_two_prime

syms th_2p_1 th_2p_2;

% length

syms l1 l2 h d;

% others

syms m1 m2 m3;

% Gravity

syms g;


%% Pipeline

l1 = 1;

l2 = 1;

start = [2; 0; 1];

goal = [1.2; 1.2; 0.5];

vmax = 1.0;

amax = 3.0;

[X, Y, Z] = Traj_samp(start, goal, vmax, amax);

[A01, A02, A03] = FWD(0, 10 * pi / 180, 0, 100, 100);

start_x = A03(1,4);
start_y = A03(2,4);

Joint_angles = INV_Traj_Jaco_prac(100, 100, [start_x;start_y;0;0;0;0], [100;100;1;0;0;0], [0;10 * pi /180;0], 100);


%% Forward Kinematics - DH Parameter

function[A01, A02, A03] = FWD(th1, th2, d, l1, l2)

    A01 = [cos(th1) -sin(th1) 0 l1*cos(th1); sin(th1) cos(th1) 0 l1*sin(th1);
        0 0 1 0; 0 0 0 1];
    
    A12 = [cos(th2) sin(th2) 0 l2*cos(th2); sin(th2) -cos(th2) 0 l2*sin(th2);
        0 0 -1 0; 0 0 0 1];
    
    A23 = [1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
    
    A02 = A01*A12;
    
    A13 = A12*A23;
    
    A03 = A01*A12*A23;
    
    px = A03(1, 4);
    
    py = A03(2, 4);
    
    pz = A03(3, 4);

end


%% Inverse Kinematics - Geometric Approach

function[Joint_angles] = INV_Geo(l1, l2, px, py, pz)

    alpha = atan2(py, px);

    R = sqrt(px^2 + py^2);

    % By Cosine 2 Law
    beta = acos((R^2 + l1^2 - l2^2) / (2*R*l1));

    th1 = alpha - beta;

    temp = (l1^2 + l2^2 - R^2) / (2*l1*l2);

    th2 = acos(-temp);

    Joint_angles = [th1; th2; 1-pz];

end


%% Inverse Kinematics - Numerical Approach (Jacobian)


function[theta] = INV_Jaco(X, Y, Z, cur_theta, l1, l2)

    dp = [X; Y; Z; 0; 0; 0];

    th1 = cur_theta(1, 1);

    th2 = cur_theta(2, 1);

    d = cur_theta(3, 1);

    J11 = -l2*sin(th1+th2) -l1*sin(th1);
        
    J12 = -l2*sin(th1+th2);
        
    J13 = 0;
        
    J21 = l2*cos(th1+th2) + l1*cos(th1);
        
    J22 = l2*cos(th1+th2);
        
    J23 = 0;
        
    J31 = 0;
        
    J32 = 0;
         
    J33 = -1;
        
    Jaco = [J11 J12 J13; J21 J22 J23; J31 J32 J33; 0 0 0; 0 0 0; 1 1 -1];
        
    psedo_inv_jaco = inv(transpose(Jaco)*Jaco)*transpose(Jaco);
        
    theta = cur_theta + psedo_inv_jaco*dp;

end


%% Dynamics

function[Torque] = DYNAMICS(th1, th2, th_p1, th_p2, th_2p_1, th_2p_2, l1, l2, h, d, m1, m2, m3, g)

    A01 = [cos(th1) -sin(th1) 0 l1*cos(th1); sin(th1) cos(th1) 0 l1*sin(th1);
            0 0 1 h; 0 0 0 1];
    
    A12 = [cos(th2) sin(th2) 0 l2*cos(th2); sin(th2) -cos(th2) 0 l2*sin(th2);
        0 0 -1 0; 0 0 0 1];
    
    A23 = [1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
    
    A02 = A01*A12;
    
    A13 = A12*A23;
    
    A03 = A01*A12*A23;

    Q1 = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];
    
    Q2 = Q1;
    
    Q3 = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
    
    
    % U
    
    U11 = [-sin(th1) -cos(th1) 0 -l1*sin(th1); cos(th1) -sin(th1) 0 l1*cos(th1);
        0 0 0 0; 0 0 0 0];
    
    U21 = [-sin(th1+th2) cos(th1+th2) 0 -l1*sin(th1)-l2*sin(th1+th2); 
        cos(th1+th2) sin(th1+th2) 0 l1*cos(th1)+l2*cos(th1+th2); 0 0 0 0; 0 0 0 0];
    
    U22 = [-sin(th1+th2) cos(th1+th2) 0 -l2*sin(th1+th2); 
        cos(th1+th2) sin(th1+th2) 0 l2*cos(th1+th2); 0 0 0 0; 0 0 0 0];
    
    U31 = [-sin(th1+th2) cos(th1+th2) 0 -l1*sin(th1)-l2*sin(th1+th2); 
        cos(th1+th2) sin(th1+th2) 0 l1*cos(th1)+l2*cos(th1+th2); 0 0 0 0; 0 0 0 0];
    
    U32 = [-sin(th1+th2) cos(th1+th2) 0 -l2*sin(th1+th2); 
        cos(th1+th2) sin(th1+th2) 0 l2*cos(th1+th2); 0 0 0 0; 0 0 0 0];
    
    U33 = [0 0 0 0; 0 0 0 0; 0 0 0 -1; 0 0 0 0];
    
    
    % J
    
    J1 = [m1*l1^2/3 0 0 -m1*l1/2; 0 0 0 0; 0 0 0 0; -m1*l1/2 0 0 m1];
    
    J2 = [m2*l2^2/3 0 0 -m2*l2/2; 0 0 0 0; 0 0 0 0; -m2*l2/2 0 0 m2];
    
    J3 = [0 0 0 0; 0 0 0 0; 0 0 m3*d^2/3 -m3*d/2; 0 0 -m3*d/2 m3];
    
    
    % Inertia Term
    
    D11 = trace(U11*J1*transpose(U11)) + trace(U21*J2*transpose(U21)) + trace(U31*J3*transpose(U31));
    
    D12 = trace(U22*J2*transpose(U21)) + trace(U32*J3*transpose(U31));
    
    D13 = trace(U33*J3*transpose(U31));
    
    D21 = trace(U21*J2*transpose(U22)) + trace(U31*J3*transpose(U32));
    
    D22 = trace(U22*J2*transpose(U22)) + trace(U32*J3*transpose(U32));
    
    D23 = trace(U33*J3*transpose(U32));
    
    D31 = trace(U31*J3*transpose(U33));
    
    D32 = trace(U32*J3*transpose(U33));
    
    D33 = trace(U33*J3*transpose(U33));
    
    
    % U - 3 parameters
    
    U111 = Q1*Q1*A01;
    
    U211 = Q1*Q1*A02;
    
    U212 = Q1*A01*Q2*A12;
    
    U221 = Q1*A01*Q2*A12;
    
    U222 = A01*Q2*Q2*A12;
    
    U311 = Q1*Q1*A03;
    
    U312 = Q1*A01*Q2*A13;
    
    U321 = Q1*A01*Q2*A13;
    
    U322 = A01*Q2*Q2*A13;
    
    
    % Colioris and Centrifugal Term
    
    % H1
    
    H111 = trace(U111*J1*transpose(U11)) + trace(U211*J2*transpose(U21)) + trace(U311*J3*transpose(U31));
    
    H112 = trace(U212*J2*transpose(U21)) + trace(U312*J3*transpose(U31));
    
    H121 = trace(U221*J2*transpose(U21)) + trace(U321*J3*transpose(U31));
    
    H122 = trace(U222*J2*transpose(U21)) + trace(U322*J3*transpose(U31));
    
    H1 = H111*th_p1*th_p1 + H112*th_p1*th_p2 + H121*th_p2*th_p1 + H122*th_p2*th_p2;
    
    % H2
    
    H211 = trace(U211*J2*transpose(U22)) + trace(U311*J3*transpose(U32));
    
    H212 = trace(U212*J2*transpose(U22)) + trace(U312*J3*transpose(U32));
    
    H221 = trace(U221*J2*transpose(U22)) + trace(U321*J3*transpose(U32));
    
    H222 = trace(U222*J2*transpose(U22)) + trace(U322*J3*transpose(U32));
    
    H2 = H211*th_p1*th_p1 + H212*th_p1*th_p2 + H221*th_p2*th_p1 + H222*th_p2*th_p2;
    
    
    % Gravity Term
    
    R11 = [-l1/2; 0; 0; 1];
    
    R22 = [-l2/2; 0; 0; 1];
    
    R33 = [0; 0; -d/2; 1];
    
    g1 = [0 0 -g 0];
    
    g2 = g1;
    
    g3 = [0 0 g 0];
    
    C1 = -m1*g1*U11*R11 -m2*g2*U21*R22 -m3*g3*U31*R33;
    
    C2 = -m2*g2*U22*R22 -m3*g3*U32*R33;
    
    C3 = -m3*g3*U33*R33;
    
    Inertia = [D11 D12 D13; D21 D22 D23; D31 D32 D33];
    
    Ang_Acc = [th_2p_1; th_2p_2; 0];
    
    Coli = [H1; H2; 0];
    
    Gravity = [C1; C2; C3];
    
    Torque = Inertia*Ang_Acc + Coli + Gravity;

end

%% Trajectory Planning - 2 1 2

function[X, Y] = Traj(start, goal, N)

    t = 7;
    Ts = 1;

    dx = goal(1, 1)-start(1, 1);

    dy = goal(2, 1)- start(2, 1);

    whole_len = sqrt(dx^2 + dy^2);

    amax = whole_len / (t*Ts - Ts^2);

    amax_dx = dx / (t*Ts - Ts^2);

    amax_dy = dy / (t*Ts - Ts^2);

    p1 = amax*Ts^2/2;
 
    p2 = amax*Ts*(t - 2*Ts);
 
    p3 = amax*Ts^2/2;

    tol = 1e-3;

    dt = t / N;

    T = 0:dt:t;
    P = ones(1, length(T));
    V = ones(1, length(T));
    A = ones(1, length(T));

    X = ones(1, length(T));
    Y = ones(1, length(T));


    for idx = 1:length(T)

        cur_t = T(idx);

        if cur_t <= Ts

            P(idx) = amax*cur_t^2/2;
            V(idx) = cur_t*amax;
            A(idx) = amax;
            X(idx) = start(1, 1) + amax_dx*cur_t^2/2;
            Y(idx) = start(2, 1) + amax_dy*cur_t^2/2;

        elseif cur_t <= t-Ts && cur_t > Ts

            cur_t = cur_t - Ts;
            P(idx) = amax*Ts^2/2 + amax*Ts*cur_t;
            V(idx) = Ts*amax;
            A(idx) = 0;
            X(idx) = start(1, 1) + amax_dx*Ts^2/2 + amax_dx*Ts*cur_t;
            Y(idx) = start(2, 1) + amax_dy*Ts^2/2 + amax_dy*Ts*cur_t;

        elseif cur_t > t- Ts && cur_t <= t

            cur_t = cur_t - (t-Ts);
            P(idx) = amax*Ts^2/2 + amax*Ts*(t-2*Ts) + cur_t*Ts*amax - amax*cur_t^2/2;
            V(idx) = (Ts-cur_t)*amax;
            A(idx) = -amax;
            X(idx) = start(1, 1) + amax_dx*Ts^2/2 + amax_dx*Ts*(t-2*Ts) + cur_t*Ts*amax_dx - amax_dx*cur_t^2/2;
            Y(idx) = start(2, 1) + amax_dy*Ts^2/2 + amax_dy*Ts*(t-2*Ts) + cur_t*Ts*amax_dy - amax_dy*cur_t^2/2;

        end

    end

    figure(4)

    subplot(3, 1, 1)
    plot(T, A, "r--")
    title("Acceleration")
    ylim([-25 25])
    hold on

    subplot(3, 1, 2)
    plot(T, V, "r--")
    title("Velocity")
    ylim([0 25])
    hold on

    subplot(3, 1, 3)
    plot(T, P, "r--")
    title("Position")
    hold on

end


%% Trajectory Planning - 2 1 2 - Sampling Time based

function[X, Y, Z] = Traj_samp(start, goal, vmax, amax)

    l1 = 1; l2 = 1;

    ts = vmax / amax;

    L = norm(goal - start);

    T = (L*amax+vmax^2)/amax*vmax;

    t = 0:0.02:T;

    disp(length(t))

    X_pre = start(1,1);
    Y_pre = start(2,1);
    Z_pre = start(3,1);

    th1 = ones(1, length(t));
    th2 = ones(1, length(t));
    d = ones(1, length(t));

    %goal_str = ['Goal : [',num2str(goal(1,1)), ",", num2str(goal(2,1)), ",",num2str(goal(3,1)), "]"];

    figure(1)
    plot3(goal(1,1), goal(2,1), goal(3,1), "r*")
    text(goal(1,1), goal(2,1), goal(3,1), 'Goal : [1.2,1.2,0.5]')
    hold on
    plot3(start(1,1), start(2,1), start(3,1),"r*")
    text(start(1,1), start(2,1), start(3,1),"Start : [2,0,1]")
    hold on

    cur_theta = zeros(3,1);

    for idx=1:length(t)

        if t(idx)<ts
            P= start+((0.5*amax*(t(idx))^2)/L)*(goal-start);

        elseif (ts<=t(idx))&&(t(idx)<(T-ts))
            P= start+((vmax*t(idx)-(vmax^2)/(2*amax))/L)*(goal-start);

        elseif t(idx)>=T-ts
            P= start+(((-amax)*(t(idx)-T)^2/2+vmax*T-vmax^2/amax)/L)*(goal-start);

        end

        X = P(1,1); Y = P(2,1); Z = P(3,1);

        cur_theta = INV_Geo(l1, l2, X, Y, Z);

        %cur_theta = INV_Jaco(X-X_pre, Y-Y_pre, Z-Z_pre, cur_theta, l1, l2);

        X_pre = X; Y_pre = Y; Z_pre = Z;

        th1(idx) = cur_theta(1,1);
        th2(idx) = cur_theta(2,1);
        d(idx) = cur_theta(3,1);

        [A01, A02, A03] = FWD(cur_theta(1, 1), cur_theta(2, 1), cur_theta(3, 1), l1, l2);

        if idx == length(t)

             figure(1)
             line = [[0;0;0], [0;0;1], [A01(1, 4); A01(2, 4); 1], [A02(1, 4); A02(2, 4); 1], [A03(1, 4);A03(2, 4);1+A03(3, 4)]];
             plot3(line(1,:),line(2,:),line(3,:),'-o')
             hold on


        else
            
            if mod(idx, 2) == 0

                figure(1)
                line = [[0;0;0], [0;0;1], [A01(1, 4); A01(2, 4); 1], [A02(1, 4); A02(2, 4); 1], [A03(1, 4);A03(2, 4);1+A03(3, 4)]];
                title('SCARA motion in Cartesian Coordinate')
                xlabel('X - axis')
                ylabel('Y - axis')
                zlabel('Z - axis')
                plot3(line(1,:),line(2,:),line(3,:),'-o')
                xlim([0 2])
                ylim([-0.5 1.5])
                zlim([0 1])
                hold on
            end
        end

    end

    disp(X_pre)
    disp(Y_pre)
    disp(Z_pre)

    Time = 1:1:length(t);

    figure(2)
    subplot(3, 1, 1)
    plot(Time, th1, 'r--', 'LineWidth', 2)
    title('joint 1 angle')
    xlim([0 length(th1)])
    hold on

    subplot(3, 1, 2)
    plot(Time, th2, 'r--', 'LineWidth', 2)
    title('joint 2 angle')
    xlim([0 length(th2)])
    hold on

    subplot(3, 1, 3)
    plot(Time, d, 'r--', 'LineWidth', 2)
    title('joint 3 dist')
    xlim([0 length(d)])
    hold on


end


%% Inverse Kinematics Using Trajectory Planning Practice (Jacobian)

function[Joint_angles] = INV_Traj_Jaco_prac(l1, l2, start, goal, cur_theta, N)

    [X, Y] = Traj(start, goal, N);

    str = 0;

    for idx = 2:N

        dp_dx = X(idx) - X(idx-1);
        dp_dy = Y(idx) - Y(idx-1);

        dp = [dp_dx; dp_dy; 0; 0; 0; 0];

        th1 = cur_theta(1,1);
        th2 = cur_theta(2,1);
        d = cur_theta(3,1);
        
        J11 = -l2*sin(th1+th2) -l1*sin(th1);
        
        J12 = -l2*sin(th1+th2);
        
        J13 = 0;
        
        J21 = l2*cos(th1+th2) + l1*cos(th1);
        
        J22 = l2*cos(th1+th2);
        
        J23 = 0;
        
        J31 = 0;
        
        J32 = 0;
         
        J33 = -1;
        
        Jaco = [J11 J12 J13; J21 J22 J23; J31 J32 J33; 0 0 0; 0 0 0; 1 1 -1];
        
        psedo_inv_jaco = inv(transpose(Jaco)*Jaco)*transpose(Jaco);
        
        next_theta = cur_theta + psedo_inv_jaco*dp;

        [A01, A02, A03] = FWD(cur_theta(1,1), cur_theta(2,1), cur_theta(3,1), l1, l2);
    
        if mod(idx, 5) == 0

            if idx == N
                disp("End")
                figure(3)
                line([0 A01(1, 4)], [0 A01(2, 4)] , 'Color', 'red')
                hold on
                line([A01(1, 4) A02(1, 4)], [A01(2, 4) A02(2, 4)])
                text(A02(1, 4)+7, A02(2, 4), "End, " + num2str(idx) + " & (" + num2str(A02(1,4)) + " , " + ...
                    num2str(A02(2, 4)) + ")", 'Color', 'blue')
                text(A02(1, 4)-50, A02(2, 4)-30, "N : " + num2str(N))
             
                hold on
                plot(A02(1, 4), A02(2, 4), "r*", 'MarkerSize', 10)

            else

                figure(3)
                title("Inverse Kinematics using Jacobian")
                str = str + 1;
                line([0 A01(1, 4)], [0 A01(2, 4)] , 'Color', 'red')
                hold on
                line([A01(1, 4) A02(1, 4)], [A01(2, 4) A02(2, 4)])
                text(A02(1, 4), A02(2, 4), num2str(idx), 'Color', 'blue')
          
                hold on

            end
    
        end

        cur_theta = next_theta;

    end

    Joint_angles = cur_theta;

end