%**************************************************************************
%
% Scenario data
% Andrea Valmorbida and Anese Giovanni
%
% This script returns the data for the estimation of position with the
% minimun squared algorithm. It gives the positions of the Ground Stations
% (GS), the initial estimation vector (position and time delay) and the
% measured pseudoranges.
%
%**************************************************************************
% receiver true position
x0_true = 4.4;
y0_true = 3.0;

switch scenario

    case 1  % no ctau
        GS_x = [9.0 0.2 3.7 -9.2 -5.5 8.4 -6.0 10.4 10.5]';   % GS x pos
        GS_y = [4.0 7.7 1.5 5.8 -2.8 -5.7 1.0 -0.8 8.0]';    % GS y pos
        toll = 3e-5;        % toll on the while loop
        k_max = 100;         % max num of iter
        ex = -0.2;          % receiver initial x pos error (sensibilità righello)
        ey = 1;         % receiver initial y pos error

        p_meas = [4.7 6.3 1.7 13.9 11.5 9.6 10.6 7.2 7.9]'; % measured pseudoranges 
        ctau = 0;        % time delay
        
        % receiver initial position
        R_x0 = x0_true + 1*ex;
        R_y0 = y0_true + 1*ey;
        
        % Initial variables values
        R_x = R_x0; R_y = R_y0; R_ctau = 0;
        Dx = [0,0]';
        k = 0;
        R_pos = [R_x,R_y];
        stop_check = 1;

    case 2  % with ctau
        GS_x = [9.0 0.2 3.7 -9.2 -5.5 8.4 -6.0 10.4 10.5]';   % GS x pos
        GS_y = [4.0 7.7 1.5 5.8 -2.8 -5.7 1.0 -0.8 8.0]';    % GS y pos
        toll = 3e-5;                    % toll on the while loop
        k_max = 100;                     % max num of iter
        ex = 0.3;                       % receiver initial x pos error
        ey = 0.5;                       % receiver initial y pos error
        
        p_meas = [4.7 6.3 1.7 13.9 11.5 9.6 10.6 7.2 7.9]'; % measured pseudoranges 
        ctau = 0.02;            % time delay
        p_meas = p_meas + ctau*0.5*unifrnd(-1,1);

        % receiver initial position
        R_x0 = x0_true + 1*ex;
        R_y0 = y0_true + 1*ey;
        R_ctau0 = ctau+ctau*unifrnd(-1,1);

        % Initial variables values
        R_x = R_x0; R_y = R_y0; R_ctau = R_ctau0;
        Dx = [0,0,0]';
        k = 0;
        R_pos = [R_x,R_y];
        stop_check = 1;

    case 3  % no ctau
        GS_x = [9.8 6.2 6.7 0.6 3.3 -2.0 -2.8]';   % GS x pos
        GS_y = [7.8 6.7 3.2 2.5 -2.7 -1.0 -6.0]';    % GS y pos
        toll = 3e-5;        % toll on the while loop, 3e-2
        k_max = 100;         % max num of iter
        ex = -0.2;          % receiver initial x pos error (sensibilità righello)
        ey = 1;         % receiver initial y pos error

        p_meas = [7.2 4.1 2.3 3.8 5.8 7.6 11.5]'; % measured pseudoranges 
        ctau = 0;        % time delay
        
        % receiver initial position
        R_x0 = x0_true + 1*ex;
        R_y0 = y0_true + 1*ey;
        
        % Initial variables values
        R_x = R_x0; R_y = R_y0; R_ctau = 0;
        Dx = [0,0]';
        k = 0;
        R_pos = [R_x,R_y];
        stop_check = 1;

    case 4  % with ctau
        GS_x = [9.8 6.2 6.7 0.6 3.3 -2.0 -2.8]';   % GS x pos
        GS_y = [7.8 6.7 3.2 2.5 -2.7 -1.0 -6.0]';    % GS y pos
        toll = 3e-5;                    % toll on the while loop
        k_max = 100;                     % max num of iter
        ex = -0.8;                       % receiver initial x pos error
        ey = 2;                       % receiver initial y pos error
        
        p_meas = [7.2 4.1 2.3 3.8 5.8 7.6 11.5]'; % measured pseudoranges 
        ctau = 0.02;            % time delay
        p_meas = p_meas + ctau*0.5*unifrnd(-1,1);

        % receiver initial position
        R_x0 = x0_true + 1*ex;
        R_y0 = y0_true + 1*ey;
        R_ctau0 = ctau+ctau*unifrnd(-1,1);
        % Initial variables values
        R_x = R_x0; R_y = R_y0; R_ctau = R_ctau0;
        Dx = [0,0,0]';
        k = 0;
        R_pos = [R_x,R_y];
        stop_check = 1;

end

% number of Ground Station (GS)
GS_n = length(GS_x);