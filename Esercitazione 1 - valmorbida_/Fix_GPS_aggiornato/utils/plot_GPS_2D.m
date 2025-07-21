%**************************************************************************
%
% Plot GPS 2D
% Andrea Valmorbida and Anese Giovanni
%
% This script plots the solutions and the covariance analysis for the GPS
% 2D Fix.
%
%**************************************************************************
%% Uncertainty ellipse
figure()

title(['GPS 2D - ', num2str(GS_n) , ' Ground Stations']); hold on

% Plot ground station
plot(GS_x,GS_y,'*','Color',[0.64 0.08 0.18],'LineWidth',2); hold on;
for j = 1:GS_n
    text(GS_x(j) + 0.5,GS_y(j),['GS ' num2str(j)],'Color',[0.64 0.08 0.18],'FontSize',16);
end

% Plot
plot(R_pos(:,1),R_pos(:,2),'o','Color',[0.1333 0.1333 0.5216],'LineWidth',1);
for j = 1:length(R_pos(:,1))
    text(R_pos(j,1),R_pos(j,2),['\bf ' num2str(j-1)],'Color',[0.1333 0.1333 0.5216],'FontSize',12);
end
plot(x0_true,y0_true,'p','Color',[0.9294    0.6941    0.1255],'LineWidth',2);
plot(xx,yy,'-','Color',[0.3922    0.8314    0.0745],'LineWidth',2);
axis equal; grid on
xlabel('X [cm]'); ylabel('Y [cm]');
ylim([min(GS_y)-2,max(GS_y)+2])
legend('','Estimated position','True Position','Uncertainty Ellipse','Location','southeast')

