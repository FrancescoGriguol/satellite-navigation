%**************************************************************************
%
% Covariance analysis
% Andrea Valmorbida and Anese Giovanni
%
% This script performs the covariance analysis on the uncertainty in
% position and the error committed during the estimation.
%
%**************************************************************************
%% Covariance on the pseudoranges
if weight == 1
    % Weighted matrix 
      Cpp = sqrt(diag(diag(eps*eps')));
      W = Cpp\eye(length(Cpp));
     
    if ctau ~= 0    % 3 unknowns
        C_xy = (H'*W*H)\eye(3);
        C_xy = C_xy(1:2,1:2);
    else            % 2 unknowns
        C_xy = (H'*W*H)\eye(2);
    end
else
    if ctau ~= 0    % 3 unknowns
        rms_2 = norm((p_meas - p_calc).^2)/(2*GS_n-3);
    else            % 2 unknowns
        rms_2 = norm((p_meas - p_calc).^2)/(2*GS_n-2);
    end

    C_xyt = rms_2 * inv(H'*H);
    C_xy = C_xyt(1:2,1:2);

end


theta = atan2(-2*C_xy(1,2),C_xy(1,1)-C_xy(2,2));
theta_deg = rad2deg(theta);
R = [sin(theta) cos(theta); -cos(theta) sin(theta)];
C_uv = R*C_xy*R';

% uncertainty ellipse
conf_val = sqrt(-2*log(1-0.997));
a = conf_val*sqrt(C_uv(1,1));
b = conf_val*sqrt(C_uv(2,2));

tt = linspace(0,2*pi,100);

uu = a*cos(tt);
vv = b*sin(tt);

uu_th =  uu*sin(theta) + vv*cos(theta);
vv_th = -uu*cos(theta) + vv*sin(theta);

xx = uu_th + R_x;
yy = vv_th + R_y;

%% Error covariance analysis
if (GS_n > 3)
    err = R_pos - [x0_true,y0_true];
    covariance = cov(err);
    [eigenvector,eigenval] = eig(covariance);
    
    smallest_eigenval = min(diag(eigenval));
    largest_eigenval = max(diag(eigenval));
    
    theta1 = atan2(eigenvector(2,2),eigenvector(1,2));
    theta1_deg = rad2deg(theta1);
    R = [sin(theta1) cos(theta1); -cos(theta1) sin(theta1)];
    
    % uncertainty ellipse
    conf_val = sqrt(-2*log(1-0.997));
    a1 = conf_val*sqrt(largest_eigenval);
    b1 = conf_val*sqrt(smallest_eigenval);
    
    tt = linspace(0,2*pi,100);
    
    uu1 = a1*cos(tt);
    vv1 = b1*sin(tt);
    
    uu_th1 =  uu1*sin(theta) + vv1*cos(theta);
    vv_th1 = -uu1*cos(theta) + vv1*sin(theta);
    
    xx1 = uu_th1 + R_x;
    yy1 = vv_th1 + R_y;
end 