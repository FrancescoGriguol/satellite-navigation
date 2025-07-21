clear 
close all

load('2022-04-12 17-20-46.bin-719743.mat');

iStart = 2500;
tStart = GPS(iStart,2)/10^6;
GPS(1:iStart,:)=[];
timeGPS = GPS(:,2)/10^6 - tStart;

IMU(IMU(1:2:end,2)/10^6<tStart,:) = []; 
XKF1(XKF1(1:2:end,2)/10^6<tStart,:) = [];

timeIMU = IMU(1:2:end,2)/10^6 - tStart;
timeKF = XKF1(1:2:end,2)/10^6 - tStart;

fIMU = 25;
fGPS = 5;
fKF = 10;
% moltiplicare la density noise density 150 micro g per la radice di dieci per
% ottenere il rumore

tGPS = 1/fGPS:1/fGPS:timeGPS(end);
tIMU = 1/fIMU:1/fIMU:timeIMU(end);
tKF = 1/fKF:1/fKF:timeKF(end);

GyrXraw = IMU(1:2:end,4);
GyrYraw = IMU(1:2:end,5);
GyrZraw = IMU(1:2:end,6);
AccXraw = IMU(1:2:end,7);
AccYraw = IMU(1:2:end,8);
AccZraw = IMU(1:2:end,9);

% Angular velocities with respect to the inertial frame coordinated in the
% body frame omwga_ib_b [rad/s]
GyrX = interp1(timeIMU,GyrXraw,tIMU);
GyrY = interp1(timeIMU,GyrYraw,tIMU);
GyrZ = interp1(timeIMU,GyrZraw,tIMU);
% Specific force vector in the body frame b [m/s^2]
AccX = interp1(timeIMU,AccXraw,tIMU);
AccY = interp1(timeIMU,AccYraw,tIMU);
AccZ = interp1(timeIMU,AccZraw,tIMU);

% latitude [deg] longitude [deg] altitude [m] of the receiver from GNSS
GPSpLat = interp1(timeGPS,GPS(:,9),tGPS);
GPSpLng = interp1(timeGPS,GPS(:,10),tGPS);
GPSpAlt = interp1(timeGPS,GPS(:,11),tGPS);

% Estimated roll, pitch and yaw angles [deg]
roll = interp1(timeKF,XKF1(1:2:end,4),tKF);
pitch = interp1(timeKF,XKF1(1:2:end,5),tKF);
yaw = interp1(timeKF,XKF1(1:2:end,6),tKF);


% Initialization

% DCM between b and g
R(:,:,1) = angle2dcm(yaw(1)*pi/180,pitch(1)*pi/180,roll(1)*pi/180);
Cg(:,:,1) = R(:,:,1)';
% Position oh the INS in the global reference frame (NED)
p(:,1) = [0 0 0]';
% Velocity of the INS in the global reference frame (NED)
pDot(:,1) = [0 0 0]';
dTIMU = 1/fIMU;
% Gravity vector in the global reference frame (NED)
gg = [0 0 9.806]';
% Transformation from the body to the global reference frame
Accg(:,1) = Cg(:,:,1)*[AccX(1) AccY(1) AccZ(1)]';

for k = 1:length(roll)-1 % k = 1 : length(GyrX)-1
    R(:,:,k) = angle2dcm(yaw(k)*pi/180,pitch(k)*pi/180,roll(k)*pi/180);
    Cg(:,:,k) = R(:,:,k)';
    Accg(:,k) = Cg(:,:,k)*[AccX(k) AccY(k) AccZ(k)]';
    pDot(:,1+k) = pDot(:,k)+(Accg(k)-gg)*dTIMU;
    p(:,1+k) = p(:,k)+pDot(:,k)*dTIMU;
end

% AVAE psi
GPSpsi = interp1(timeGPS,GPS(:,13),tGPS);
Accpsi = zeros(3,length(GPSpsi));
for i = 1:length(GPSpsi)
    Accpsi(1,i) = Accg(1,i)*cos(GPSpsi(i)*pi/180)+Accg(2,i)*sin(GPSpsi(i)*pi/180);
    Accpsi(2,i) = Accg(2,i)*cos(GPSpsi(i)*pi/180)-Accg(1,i)*sin(GPSpsi(i)*pi/180);
    Accpsi(3,i) = Accg(3,i)-9.806;
end

% AVAE theta
AccXpsi = Accpsi(1,:)';
AccYpsi = Accpsi(2,:)';
AccZpsi = Accpsi(3,:)';
INStheta = zeros(1,length(GPSpsi));
for i = 1:length(GPSpsi)
    num = AccXpsi(i)*AccZpsi(i)+AccXpsi(i)*sqrt(AccXpsi(i)^2+AccZpsi(i)^2-AccXpsi(i)^2);
    den = AccZpsi(i)^2-AccXpsi(i)^2;
    INStheta(1,i) = atan2(num,den);
end

% AVAE roll



figure()
geoplot(GPS(:,9),GPS(:,10),'.')
geobasemap satellite

figure
plot(tKF,roll)
hold on
grid on
plot(tKF,pitch)
plot(tKF,yaw)
legend('Roll','Pitch','Yaw')
title('Extended Kalman Filter euler angles')
ylabel('\phi,\theta,\psi [deg]')
xlabel('time [s]')

figure
plot(tIMU,rollINS*180/pi)
hold on
grid on
plot(tIMU,pitchINS*180/pi)
plot(tIMU,yawINS*180/pi)

figure
plot(tKF,yaw)
hold on
grid on
plot(tIMU,yawINS*180/pi)

figure
grid on 
hold on
plot(tIMU,AccX)
plot(tIMU,AccY)
plot(tIMU,AccZ)

figure
grid on 
hold on
plot(tIMU(1:end-1),Acc_e(1,:))
plot(tIMU(1:end-1),Acc_e(2,:))
plot(tIMU(1:end-1),Acc_e(3,:))