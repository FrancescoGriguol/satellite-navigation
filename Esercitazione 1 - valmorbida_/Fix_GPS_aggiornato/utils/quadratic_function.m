% Quadratic function
% Chi squared
% x,y meshgrid
[x,y]=meshgrid(linspace(-max(abs([GS_x', GS_y'])),max(abs([GS_x', GS_y'])),100));
x=x+R_x; y=y+R_y;

% compute the chi squared value for each pair (x,y)
z = sqrt((x-R_x).^2+(y-R_y).^2)+ctau;
z_GS = sqrt((GS_x-R_x).^2+(GS_y-R_y).^2);
z_i = sqrt((R_pos(:,1)-R_x).^2+(R_pos(:,2)-R_y).^2); % i-th chi squared value for the i-th iteration

figure(1)
surf(x,y,z,"FaceAlpha",0.5,'LineStyle','none')
hold on
contour3(x,y,z,0.5:0.5:10,'-','LineWidth',1)
scatter3(GS_x,GS_y,z_GS,100,"rs")

% plot the points for the chi squared values of the iterations
scatter3(R_pos(:,1),R_pos(:,2), z_i, 100, 'r+','LineWidth',2)

hold off

xlabel('x [cm]')
ylabel('y [cm]')
zlabel('\chi^2 [cm^2]')

% Sigma squared
delta_x = linspace(-20,20,100);
delta_y = linspace(-20,20,100);

% preallocation for the sigma squared matrix 
z_sigma = zeros(length(delta_x));

if ctau == 0    % 2 unknwons
    for i = 1:length(delta_x)
        for j = 1:length(delta_y)
            z_sigma(i,j) = (eps-H*[delta_x(i);delta_y(j)])'*(eps-H*[delta_x(i);delta_y(j)]);
        end
    end
    
else    
    for i = 1 : length(delta_x)
        for j = 1 : length(delta_y)
            z_sigma(i,j) = (eps-H*[delta_x(i);delta_y(j);0])'*(eps-H*[delta_x(i);delta_y(j);0]);
        end
    end
end

% sigma squared value for the k-th dx vector
z_i(k) = (eps-H*Dx)'*(eps-H*Dx);
Dx_tot = [Dx_tot, Dx];

figure(2)
surf(delta_x,delta_y,z_sigma,"FaceAlpha",0.2,'LineStyle','none')
hold on
contour3(delta_x,delta_y,z_sigma,0:100:5000,'-','LineWidth',1)
scatter3(Dx(1,end),Dx(2,end),z_i(1:end),100,'filled',"o")
view([-70, 45])
hold off
xlabel('Dx [cm]')
ylabel('Dy [cm]')
zlabel('\Sigma^2 [cm^2]')
