function [distance, yy_bond_sc_correlation, yy_prime_bond_sc_correlation] = read_hybrid_space_SC(FileNamePostfix, Ly)

directory_to_hubbard = '/home/haoxinwang/GitHub/Hubbard/';
%%===== Preprocess =====%%
SCData = jsondecode(fileread([directory_to_hubbard, '/data/onsitepair00',FileNamePostfix]));
group_data_num = numel(SCData);
distance = zeros(1,group_data_num);
for i = 1:group_data_num
    site1 = SCData{i}{1}(1); site2 = SCData{i}{1}(2);
    distance(i) = ((site2)-(site1))/Ly;
end
zeromomentum_yy_bond_sc_correlation = zeros(1,group_data_num); % zero momentum in y-direction
yy_bond_sc_correlation = zeros(1,group_data_num); % \Phi_{yy}(r)
yy_prime_bond_sc_correlation = zeros(1,group_data_num); % \Phi_{yy}(r)
%%===== End of preprocess =====%%


%%===== Group 1: 2-point correlation data in hybrid space =====%%
for k1 = 0:Ly-1
    for k2 = 0:Ly-1
        if( mod(2*(k2-k1)+Ly, Ly)==0)
            SCData = jsondecode(fileread([directory_to_hubbard,'data/onsitepair',num2str(k1),num2str(k2),FileNamePostfix]));
            onsite_pair_raw_data = zeros(1,group_data_num);
            for i = 1:group_data_num
                %                 site1 = SCData{i}{1}(1); site2 = SCData{i}{1}(2);
                %                 distance(i) = ((site2-k2)-(site1-k1))/Ly;
                onsite_pair_raw_data(i) = SCData{i}{2};
            end
            K1 = k1 * 2*pi/Ly; K2 = k2 * 2 * pi /Ly;
            if (mod(2*k1,Ly)==0)
                zeromomentum_yy_bond_sc_correlation = zeromomentum_yy_bond_sc_correlation + 2*cos(K1)*cos(K2)*onsite_pair_raw_data;
            end
            yy_bond_sc_correlation = yy_bond_sc_correlation + 4 * cos(K1-K2) * onsite_pair_raw_data;
            yy_prime_bond_sc_correlation = yy_prime_bond_sc_correlation + 4 * cos(K1-K2) * cos(2*K2) * onsite_pair_raw_data;
            %loglog(distance, abs(onsite_pair_raw_data),'-x');hold on;
        end
    end
end
%%===== End of Group 1: 2-point correlation in hybrid space =====%%

% for 4-leg the following two parts are actually no contribute
%%===== Group 2: 3-point correlation (A) data in hybrid space =====%%
for k1 = 0:Ly-1
    for k2 = 0:Ly-2
        for k3 = k2+1:Ly-1
            if ( mod ( 2*k1 - k2 - k3 + 2*Ly, Ly)==0)
                SCData_a = jsondecode(fileread([directory_to_hubbard,'data/scsA',num2str(k1),num2str(k2),num2str(k3),'a',FileNamePostfix]));
                SCData_b = jsondecode(fileread([directory_to_hubbard,'data/scsA',num2str(k1),num2str(k2),num2str(k3),'b',FileNamePostfix]));
                %                 fprintf("k1 = %d, k2 = %d, k3 = %d\n",k1,k2,k3);
                pair3A_raw_data = zeros(1, group_data_num);
                for i = 1:group_data_num
                    pair3A_raw_data(i) = -SCData_a{i}{2} -SCData_b{i}{2};
                end
                K1 = k1 * 2*pi/Ly; K2 = k2 * 2 * pi /Ly; K3 = k3*2*pi/Ly;
                if (mod(2*k1,Ly)==0)
                    zeromomentum_yy_bond_sc_correlation = zeromomentum_yy_bond_sc_correlation + 2*cos(K1)*cos(K2)*pair3A_raw_data;
                end
                yy_bond_sc_correlation = yy_bond_sc_correlation + 2 * (cos(K1-K2) + cos(K1-K3)) * pair3A_raw_data; %actually these two lines can use the same form with previous 2-point function.
                yy_prime_bond_sc_correlation = yy_prime_bond_sc_correlation + 2 * (cos(K1-K2) + cos(K1-K3)) * cos(K2+K3) * pair3A_raw_data;
            end
        end
    end
end
%%===== End of Group 2: 3-point correlation (A) data in hybrid space =====%%


%%===== Group 3: 3-point correlation (B) data in hybrid space =====%%
for k1 = 0:Ly-1
    for k2 = 0:Ly-2
        for k3 = k2+1:Ly-1
            if ( mod ( 2*k1 - k2 - k3 + 2*Ly, Ly)==0)
                SCData_a = jsondecode(fileread([directory_to_hubbard,'data/scsB',num2str(k1),num2str(k2),num2str(k3),'a',FileNamePostfix]));
                SCData_b = jsondecode(fileread([directory_to_hubbard,'data/scsB',num2str(k1),num2str(k2),num2str(k3),'b',FileNamePostfix]));
                
                pair3B_raw_data = zeros(1, group_data_num);
                for i = 1:group_data_num
                    pair3B_raw_data(i) = -SCData_a{i}{2} -SCData_b{i}{2};
                end
                K1 = k1 * 2*pi/Ly; K2 = k2 * 2 * pi /Ly; K3 = k3*2*pi/Ly;
                if (mod(2*k1,Ly)==0)
                    zeromomentum_yy_bond_sc_correlation = zeromomentum_yy_bond_sc_correlation + 2*cos(K1)*cos(K2)*pair3B_raw_data;
                end
                yy_bond_sc_correlation = yy_bond_sc_correlation + 2 * (cos(K1-K2) + cos(K1-K3)) * pair3B_raw_data;
                yy_prime_bond_sc_correlation = yy_prime_bond_sc_correlation + 2 * (cos(K1-K2) + cos(K1-K3)) * cos(K2+K3) * pair3B_raw_data;
            end
        end
    end
end
%%===== End of Group 3: 3-point correlation (B) data in hybrid space =====%%

%%===== Group 4: 4-point correlation data in hybrid space =====%%
for k1 = 0:Ly-2
    for k2 = k1+1:Ly-1
        for k3 = 0:Ly-2
            for k4 = k3+1:Ly-1
                if mod(k1+k2-k3-k4 + 2* Ly,Ly)==0
                    SCData_a = jsondecode(fileread([directory_to_hubbard,'data/scs',num2str(k1),num2str(k2),num2str(k3),num2str(k4),'a',FileNamePostfix]));
                    SCData_b = jsondecode(fileread([directory_to_hubbard,'data/scs',num2str(k1),num2str(k2),num2str(k3),num2str(k4),'b',FileNamePostfix]));
                    SCData_c = jsondecode(fileread([directory_to_hubbard,'data/scs',num2str(k1),num2str(k2),num2str(k3),num2str(k4),'c',FileNamePostfix]));
                    SCData_d = jsondecode(fileread([directory_to_hubbard,'data/scs',num2str(k1),num2str(k2),num2str(k3),num2str(k4),'d',FileNamePostfix]));
                    pair4_raw_data = zeros(1, group_data_num);
                    for i = 1:group_data_num
                        pair4_raw_data(i) = SCData_a{i}{2} + SCData_b{i}{2}+ SCData_c{i}{2}+ SCData_d{i}{2};
                    end
                    K1 = k1 * 2*pi/Ly; K2 = k2 * 2 * pi /Ly; K3 = k3*2*pi/Ly;  K4 = k4*2*pi/Ly;
                    if(mod(k1+k2, Ly) ==0)
                        zeromomentum_yy_bond_sc_correlation = zeromomentum_yy_bond_sc_correlation + 2*cos(K1)*cos(K3)*pair4_raw_data; %in 4-leg, useless
                    end
                    yy_bond_sc_correlation = yy_bond_sc_correlation + 2 * (cos(K1-K3) + cos(K1-K4)) * pair4_raw_data; 
                    yy_prime_bond_sc_correlation = yy_prime_bond_sc_correlation + 2 * (cos(K1-K3) + cos(K1-K4)) * cos(K3+K4) * pair4_raw_data; %in 4-leg, useless
                end
            end
        end
    end
end
%%===== End of Group 4: 4-point correlation data in hybrid space =====%%




end