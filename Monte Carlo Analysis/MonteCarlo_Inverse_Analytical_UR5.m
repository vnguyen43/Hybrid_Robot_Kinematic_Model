clear all
close all
clc

cd('C:\Users\vtn2\Documents\MATLAB - Local\z_Coding Experiments');

%Transformation: Robot = r2s*TRIO
r2s = csvread('Robot2Sensor_UR5.csv');

data = csvread('Data_UR5_0250.csv');
data( ~any(data,2), : ) = []; %Remove 0's

q1 = data(:,1); q2 = data(:,2); q3 = data(:,3); q4 = data(:,4); q5 = data(:,5); q6 = data(:,6);
q = [q1, q2, q3, q4, q5, q6];
px = data(:,7); py = data(:,8); pz = data(:,9); rx = deg2rad(data(:,10)); ry = deg2rad(data(:,11)); rz = deg2rad(data(:,12));

q1_RMSE = 0.0720000000000000; %Root mean squared error from inverse evaluation (degrees)
numSim = 2; %Number of simulations for each MC

%Factory Calibrated Correction Factor
delta_a = [0.00022665,0.340369,0.00010277,-0.0003363,0.00045491,0]; %mm?
delta_d = [0.0001078,380.7757,-381.459,0.684885,-2.2636e-5,4.01149e-5];%mm?
delta_p = [-0.0002813,0.0010939,0.0079698,-0.00082227,0.000452722,0]; %radians?

%Dh parameters
a0 = 0;     a1 = 0+delta_a(1);             a2 = -0.425*1000+delta_a(2);     a3 = -0.39225*1000+delta_a(3);  a4 = 0+delta_a(4);                 a5 = 0+delta_a(5);             a6 = 0+delta_a(6);% a (unit: mm)
d0 = 0;     d1 = 0.08916*1000+delta_d(1);   d2 = 0+delta_d(2);               d3 = 0+delta_d(3);             d4 = 0.10915*1000+delta_d(4);     d5 = 0.09465*1000+delta_d(5);   d6 = 0.0823*1000+delta_d(6);% d (unit: mm)
p0 = 0;     p1 = +pi/2+delta_p(1);         p2 = 0+delta_p(2);               p3 = 0+delta_p(3);             p4 = +pi/2+delta_p(4);             p5 = -pi/2+delta_p(5);         p6 = 0; %alpha


for g = 1:1:length(q1)
    g
    
    k = 1;
    
    trio = eye(4);
    trio(1:3,4) = [px(g);py(g);pz(g)];
    trio(1:3,1:3) = (SpinCalc('EA321toDCM',rad2deg([rz(g),ry(g),rx(g)]),0))';
    H06 = r2s*trio;
    
    H60 = H06^-1;
    
    q1_a(g,1) = deg2rad(q1(g));
    q1_mc(k,:) = (q1_a(g,1) + randn(numSim,1)*deg2rad(q1_RMSE))';
    
    q5_mc(k,:) = +(acos((H06(1,4)*sin(q1_mc(k,:))-H06(2,4)*cos(q1_mc(k,:))-d4)/d6));
    
    q6_temp = [];
    for h = 1:1:numSim
        temp = ((atan2((-H60(2,1)*sin(q1_mc(k,:))+H60(2,2)*cos(q1_mc(k,:)))/sin(q5_mc(k,h)),H60(1,1)*sin(q1_mc(k,:))-H60(1,2)*cos(q1_mc(k,:))/sin(q5_mc(k,h)))));
        q6_temp = [q6_temp, temp];
    end
    q6_mc(k,:) = q6_temp;
    
    q3_temp = [];
    for a=1:1:length(q1_mc(k,:))
        for b=1:1:length(q5_mc(k,:))
            for c=1:1:length(q6_mc(k,:))
                H01 = [
                    cos(q1_mc(k,a)),            -sin(q1_mc(k,a)),           0,          a0; ...
                    sin(q1_mc(k,a))*cos(p0),    cos(q1_mc(k,a))*cos(p0),    -sin(p0),    -sin(p0)*d1; ...
                    sin(q1_mc(k,a))*sin(p0),    cos(q1_mc(k,a))*sin(p0),    cos(p0),    +cos(p0)*d1; ...
                    0,                  0,                   0,           1
                    ];
                
                H45 = [
                    cos(q5_mc(k,b)),            -sin(q5_mc(k,b)),           0,          a4; ...
                    sin(q5_mc(k,b))*cos(p4),    +cos(q5_mc(k,b))*cos(p4),    -sin(p4),    -sin(p4)*d5; ...
                    sin(q5_mc(k,b))*sin(p4),    +cos(q5_mc(k,b))*sin(p4),    +cos(p4),    +cos(p4)*d5; ...
                    0,                  0,                   0,           1
                    ];
                
                H56 = [
                    cos(q6_mc(k,c)),            -sin(q6_mc(k,c)),           0,          a5; ...
                    sin(q6_mc(k,c))*cos(p5),    +cos(q6_mc(k,c))*cos(p5),    -sin(p5),    -sin(p5)*d6; ...
                    sin(q6_mc(k,c))*sin(p5),    +cos(q6_mc(k,c))*sin(p5),    +cos(p5),    +cos(p5)*d6; ...
                    0,                  0,                   0,           1
                    ];
                
                H14 = H01^-1*H06*(H45*H56)^-1;
                norm_P14 = sqrt(H14(1,4)^2+H14(2,4)^2+H14(3,4)^2);
                temp = +((acos((norm_P14^2-a2^2-a3^2)/(2*a2*a3))));
                q3_temp = [q3_temp, temp];
            end
        end
    end
    q3_mc(k,:) = q3_temp;
    
    q2_temp = [];
    for a=1:1:length(q1_mc(k,:))
        for b=1:1:length(q5_mc(k,:))
            for c=1:1:length(q6_mc(k,:))
                for d=1:1:length(q3_mc(k,:))
                    H01 = [
                        cos(q1_mc(k,a)),            -sin(q1_mc(k,a)),           0,          a0; ...
                        sin(q1_mc(k,a))*cos(p0),    cos(q1_mc(k,a))*cos(p0),    -sin(p0),    -sin(p0)*d1; ...
                        sin(q1_mc(k,a))*sin(p0),    cos(q1_mc(k,a))*sin(p0),    cos(p0),    +cos(p0)*d1; ...
                        0,                  0,                   0,           1
                        ];
                    
                    H45 = [
                        cos(q5_mc(k,b)),            -sin(q5_mc(k,b)),           0,          a4; ...
                        sin(q5_mc(k,b))*cos(p4),    +cos(q5_mc(k,b))*cos(p4),    -sin(p4),    -sin(p4)*d5; ...
                        sin(q5_mc(k,b))*sin(p4),    +cos(q5_mc(k,b))*sin(p4),    +cos(p4),    +cos(p4)*d5; ...
                        0,                  0,                   0,           1
                        ];
                    
                    H56 = [
                        cos(q6_mc(k,c)),            -sin(q6_mc(k,c)),           0,          a5; ...
                        sin(q6_mc(k,c))*cos(p5),    +cos(q6_mc(k,c))*cos(p5),    -sin(p5),    -sin(p5)*d6; ...
                        sin(q6_mc(k,c))*sin(p5),    +cos(q6_mc(k,c))*sin(p5),    +cos(p5),    +cos(p5)*d6; ...
                        0,                  0,                   0,           1
                        ];
                    
                    H14 = H01^-1*H06*(H45*H56)^-1;
                    norm_P14 = sqrt(H14(1,4)^2+H14(2,4)^2+H14(3,4)^2);
                    
                    temp = (atan2(-H14(3,4),-H14(1,4))-asin(-a3*sin(q3_mc(k,d))/norm_P14));
                    q2_temp = [q2_temp, temp];
                end
            end
        end
    end
    
    q2_mc(k,:) = q2_temp;
    
    q4_temp = [];
    for a=1:1:length(q1_mc(k,:))
        for b=1:1:length(q5_mc(k,:))
            for c=1:1:length(q6_mc(k,:))
                for d=1:1:length(q3_mc(k,:))
                    for e=1:1:length(q2_mc(k,:))
                        H01 = [
                            cos(q1_mc(k,a)),            -sin(q1_mc(k,a)),           0,          a0; ...
                            sin(q1_mc(k,a))*cos(p0),    cos(q1_mc(k,a))*cos(p0),    -sin(p0),    -sin(p0)*d1; ...
                            sin(q1_mc(k,a))*sin(p0),    cos(q1_mc(k,a))*sin(p0),    cos(p0),    +cos(p0)*d1; ...
                            0,                  0,                   0,           1
                            ];
                        
                        H45 = [
                            cos(q5_mc(k,b)),            -sin(q5_mc(k,b)),           0,          a4; ...
                            sin(q5_mc(k,b))*cos(p4),    +cos(q5_mc(k,b))*cos(p4),    -sin(p4),    -sin(p4)*d5; ...
                            sin(q5_mc(k,b))*sin(p4),    +cos(q5_mc(k,b))*sin(p4),    +cos(p4),    +cos(p4)*d5; ...
                            0,                  0,                   0,           1
                            ];
                        
                        H56 = [
                            cos(q6_mc(k,c)),            -sin(q6_mc(k,c)),           0,          a5; ...
                            sin(q6_mc(k,c))*cos(p5),    +cos(q6_mc(k,c))*cos(p5),    -sin(p5),    -sin(p5)*d6; ...
                            sin(q6_mc(k,c))*sin(p5),    +cos(q6_mc(k,c))*sin(p5),    +cos(p5),    +cos(p5)*d6; ...
                            0,                  0,                   0,           1
                            ];
                        
                        H14 = H01^-1*H06*(H45*H56)^-1;
                        H12 = [
                            cos(q2_mc(k,e)),            -sin(q2_mc(k,e)),           0,          a1; ...
                            sin(q2_mc(k,e))*cos(p1),    +cos(q2_mc(k,e))*cos(p1),    -sin(p1),    -sin(p1)*d2; ...
                            sin(q2_mc(k,e))*sin(p1),    +cos(q2_mc(k,e))*sin(p1),    +cos(p1),    +cos(p1)*d2; ...
                            0,                  0,                   0,           1
                            ];
                        
                        H23 = [
                            cos(q3_mc(k,d)),            -sin(q3_mc(k,d)),           0,          a2; ...
                            sin(q3_mc(k,d))*cos(p2),    +cos(q3_mc(k,d))*cos(p2),    -sin(p2),    -sin(p2)*d3; ...
                            sin(q3_mc(k,d))*sin(p2),    +cos(q3_mc(k,d))*sin(p2),    +cos(p2),    +cos(p2)*d3; ...
                            0,                  0,                   0,           1
                            ];
                        
                        H34 = (H12*H23)^-1*H14;
                        temp = (atan2(H34(2,1),H34(1,1)));
                        q4_temp = [q4_temp, temp];
                    end
                end
            end
        end
    end
%     q4_mc(k,:) = q4_temp;

    q_std(g,:) = rad2deg([std(q1_mc(k,:)),std(q2_mc(k,:)),std(q3_mc(k,:)),std(q4_temp),std(q5_mc(k,:)),std(q6_mc(k,:))]);
    
end