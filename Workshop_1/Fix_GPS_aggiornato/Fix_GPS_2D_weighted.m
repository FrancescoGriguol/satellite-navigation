%% GPS 2D
% code developed by:
% Andrea Valmorbida and Giovanni Anese

% clear workspace
%clear; clc; close;

addpath("utils/")

% set Default Axes FontSize to 18
set(0,'DefaultAxesFontSize',24);

% set the scenario
% 1: no ctau and no aligned GS
% 2: ctau and no aligned GS
% 3: no ctau and aligned GS
% 4: ctau and aligned GS

scenario = 2; 
scenario_data;

% select if weighted (1) or not (0)
weight = 0;

%%


% Ground Station ID
GS_num = (1:GS_n)';

            
fprintf('%6s %14s %14s %14s %14s \n', 'iter', 'dx [cm]', 'dy [cm]', 'x [cm]', 'y [cm]')
fprintf('%6i %14.6f %14.6f %14.6f %14.6f \n',[0,0,0,R_x0,R_y0])
%% main while loop

% memalloc
d_calc = zeros(GS_n,1);
p_calc = zeros(GS_n,1);

if ctau ~= 0  
    H = zeros(GS_n,3);
else 
    H = zeros(GS_n,2);
end
Dx_tot = [];

while ((stop_check>toll) && (k<k_max))
    
    % d_calc
    for i=1:GS_n
        d_calc(i,1) = sqrt((GS_x(i)-R_x)^2+(GS_y(i)-R_y)^2);
        p_calc(i,1) = d_calc(i,1) + R_ctau;
    end
    
    % eps
    eps = p_meas - p_calc;
    
    % H matrix
    for i=1:GS_n
        H(i,1) = -(GS_x(i)-R_x)/d_calc(i);
        H(i,2) = -(GS_y(i)-R_y)/d_calc(i);
        if ctau ~= 0
            H(i,3) = 1;
        end
    end

    k = k + 1;
    
    % Weight matrix
    if weight == 1
        Cpp = sqrt(diag(diag(eps*eps')));
        W = Cpp\eye(length(Cpp));
    else
        W = eye(GS_n);
    end

    if (GS_n == 3 && ctau ~= 0)
        Dx_new = (W*H)\eps;
    else
        AA = H'*W*H;
        bb = H'*W*eps;
        
        Dx_new = AA\bb; 
    end
    
    % stop condition
    stop_check = norm(Dx_new(1:2));
    
    % receiver position
    Dx = Dx_new;
    R_x = R_x + Dx(1);
    R_y = R_y + Dx(2);
    if length(Dx) == 3
        R_ctau = R_ctau + Dx(3);
    end
    
    R_pos = [R_pos;R_x,R_y];
    
    fprintf('%6i %14.6e %14.6e %14.6f %14.6f \n',[k,Dx(1),Dx(2),R_x,R_y])

    % Quadratic function
    quadratic_function;

end
fprintf('\n')
fprintf('Receiver true position      = [%.3f, %.3f] cm\n',x0_true,y0_true);
fprintf('Receiver estimated position = [%.3f, %.3f] cm\n',R_x,R_y);
fprintf('Estimation error            = [%.3f, %.3f] cm\n',R_x-x0_true,R_y-y0_true);


%% Covariance analysis

covariance_analysis_weighted;

%% Plot

plot_GPS_2D;