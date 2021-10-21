clear all
close all
clc

cd('C:\Users\vtn2\Documents\MATLAB - Local\z_Coding Experiments');

%Transformation: Robot = r2s*TRIO
r2s = csvread('Robot2Sensor_UR5.csv');

data = csvread('Data_UR5_0250.csv');
data( ~any(data,2), : ) = []; %Remove 0's

q_deg = data(:,1:6);
q = deg2rad(q_deg);

q1_rad = q(:,1); q2_rad = q(:,2); q3_rad = q(:,3); q4_rad = q(:,4); q5_rad = q(:,5); q6_rad = q(:,6);
q1_deg = q_deg(:,1); q2_deg = q_deg(:,2); q3_deg = q_deg(:,3); q4_deg = q_deg(:,4); q5_deg = q_deg(:,5); q6_deg = q_deg(:,6);

px = data(:,7); py = data(:,8); pz = data(:,9); rx = deg2rad(data(:,10)); ry = deg2rad(data(:,11)); rz = deg2rad(data(:,12));

load('Hybrid_Mdl_UR5');

for g = 1:1:length(q1_deg)
    g
% for g = 1:1:1
    
    trio = eye(4);
    trio(1:3,4) = [px(g);py(g);pz(g)];
    trio(1:3,1:3) = (SpinCalc('EA321toDCM',rad2deg([rz(g),ry(g),rx(g)]),0))';
    H06 = r2s*trio;
    
    H60 = H06^-1;
    
    H_input(g,:) = [H06(1,1),H06(1,2),H06(1,3),H06(1,4),H06(2,1),H06(2,2),H06(2,3),H06(2,4),H06(3,1),H06(3,2),H06(3,3),H06(3,4)];
    
    
    q1_RMSE = 0.05500000; %Root mean squared error from inverse evaluation (degrees)
    numSim = 2; %Number of simulations for each MC
    
    q1_a = deg2rad(q1_deg(g));
    q1_mc(1,:) = (q1_a + randn(numSim,1)*deg2rad(q1_RMSE))';
    
    for a = 1:1:length(q1_mc)
        q5_mc(1,a) = deg2rad(predict(h_q5,[H_input(g,:),q1_mc(a)]));
    end
    
    ind = 1; q6_mc = zeros(1,length(q1_mc)*length(q5_mc));
    for a =1:1:length(q1_mc)
        for b= 1:1:length(q5_mc)
        q6_mc(1,ind) = deg2rad(predict(h_q6,[H_input(g,:),q1_mc(a),q5_mc(b)]));
        ind = ind + 1;
        end
    end
    
    ind = 1; q3_mc = zeros(1,length(q1_mc)*length(q5_mc)*length(q6_mc));
    for a =1:1:length(q1_mc)
        for b= 1:1:length(q5_mc)
            for c= 1:1:length(q6_mc)
            q3_mc(1,ind) = deg2rad(predict(h_q3,[H_input(g,:),q1_mc(a),q5_mc(b),q6_mc(c)]));
            ind = ind + 1;
            end
        end
    end
    
    ind = 1; q2_mc = zeros(1,length(q1_mc)*length(q5_mc)*length(q6_mc)*length(q3_mc));
    for a =1:1:length(q1_mc)
        for b= 1:1:length(q5_mc)
            for c= 1:1:length(q6_mc)
                for d= 1:1:length(q3_mc)
                q2_mc(1,ind) = deg2rad(predict(h_q2,[H_input(g,:),q1_mc(a),q5_mc(b),q6_mc(c),q3_mc(d)]));
                ind = ind + 1;
                end
            end
        end
    end
    
    ind = 1; q4_mc = zeros(1,length(q1_mc)*length(q5_mc)*length(q6_mc)*length(q3_mc)*length(q2_mc));
    for a =1:1:length(q1_mc)
        for b= 1:1:length(q5_mc)
            for c= 1:1:length(q6_mc)
                for d= 1:1:length(q3_mc)
                    for e =1:1:length(q2_mc)
                    q4_mc(1,ind) = deg2rad(predict(h_q4,[H_input(g,:),q1_mc(a),q5_mc(b),q6_mc(c),q3_mc(d),q2_mc(e)]));
                    ind = ind + 1;
                    end
                end
            end
        end
    end
    
    q_std(g,:) = rad2deg([std(q1_mc),std(q2_mc),std(q3_mc),std(q4_mc),std(q5_mc),std(q6_mc)]);

    
end

save('MC_Hybrid_UR5')
