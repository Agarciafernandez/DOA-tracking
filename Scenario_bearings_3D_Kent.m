randn('seed',9)
rand('seed',9)


Nsteps=100; %Time steps in the simulation
Nmc=1000; %Number of Monte Carlo runs
Nmc_trajectory=20; %Number of Monte Carlo runs for a fixed trajectories
N_trajectories=Nmc/Nmc_trajectory; %Number of trajectories


%Prior of trajectory
x0=[-100;5;0;5;50;4];
P_ini=diag([5^2,1,150^2,1,4^2,0.1^2]);




T=0.5;
sigmaU=0.5;

F_single=[1 T;0 1];
Q_single=sigmaU^2*[T^3/3 T^2/2;T^2/2 T];

Q=kron(eye(3),Q_single);
F=kron(eye(3),F_single);




%Measurement model

%Sensor positions
x_s_1=[100;0;0];
x_s_2=[-250;0;0];

%Sensor noise
sigma2_bear=(3*pi/180)^2;
sigma2_el=(1*pi/180)^2;


if(sigma2_bear<sigma2_el)
    display('Error in choice of sigma bearings/elevation')
end
kappa=1/2*(1/sigma2_el+1/sigma2_bear);
beta=1/4*(1/sigma2_el-1/sigma2_bear);



measurement_kent=1; %If it is 1, measurements are generated from Kent distribution, otherwise from Gaussian


chol_ini=chol(P_ini)';
chol_Q=chol(Q)';



X_multi_series=zeros(N_trajectories,6,Nsteps);


for i=1:N_trajectories
    X_multi_i=zeros(6,Nsteps);
    xk=x0+chol_ini*randn(6,1);
    X_multi_i(:,1)=xk;
    for k=1:Nsteps-1
        xk_pred=F*xk+chol_Q*randn(6,1);
        X_multi_i(:,k+1)=xk_pred;
        xk=xk_pred;
    end
    X_multi_series(i,:,:)=X_multi_i;
end

% %Plot scenario
figure(1)
clf
for i=1:N_trajectories
    X_multi_i=squeeze(X_multi_series(i,:,:));
    plot3(X_multi_i(1,:),X_multi_i(3,:),X_multi_i(5,:),'b')
    hold on
end
plot3(x_s_1(1),x_s_1(2),x_s_1(3),'xr','Linewidth',3,'MarkerSize',10)
plot3(x_s_2(1),x_s_2(2),x_s_2(3),'xr','Linewidth',3,'MarkerSize',10)
grid on
axis equal
hold off
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
% 
% %Plot angle-bearing trajectories
% figure(2)
% clf
% for i=1:N_trajectories
%     X_multi_i=squeeze(X_multi_series(i,:,:));
%     dist_sensor=sqrt((X_multi_i(1,:)-x_s_1(1)).^2+(X_multi_i(3,:)-x_s_1(2)).^2+(X_multi_i(5,:)-x_s_1(3)).^2);
%     plot(180/pi*atan2(X_multi_i(3,:)-x_s_1(2),X_multi_i(1,:)-x_s_1(1)),180/pi*asin((X_multi_i(5,:)-x_s_1(3))./dist_sensor),'-b','Linewidth',1.3)
%     hold on
%     dist_sensor=sqrt((X_multi_i(1,:)-x_s_2(1)).^2+(X_multi_i(3,:)-x_s_2(2)).^2+(X_multi_i(5,:)-x_s_2(3)).^2);
%     plot(180/pi*atan2(X_multi_i(3,:)-x_s_2(2),X_multi_i(1,:)-x_s_2(1)),180/pi*asin((X_multi_i(5,:)-x_s_2(3))./dist_sensor),'--r','Linewidth',1.3)
%     
% end
% grid on
% %axis equal
% axis(180/pi*[-pi pi, 0 pi/2])
% hold off
% xlabel('Azimuth (º)')
% ylabel('Elevation (º)')
% legend('Sensor 1','Sensor 2')

%Plot distance to sensors

% figure(3)
% clf
% average_dist_sensor1=zeros(1,Nsteps);
% average_dist_sensor2=zeros(1,Nsteps);
% 
% for i=1:N_trajectories
%     X_multi_i=squeeze(X_multi_series(i,:,:));
%     dist_sensor=sqrt((X_multi_i(1,:)-x_s_1(1)).^2+(X_multi_i(3,:)-x_s_1(2)).^2+(X_multi_i(5,:)-x_s_1(3)).^2);
%     average_dist_sensor1=average_dist_sensor1+dist_sensor;
%     plot(1:Nsteps,dist_sensor,'-b','Linewidth',1.3)
%     hold on
%     dist_sensor=sqrt((X_multi_i(1,:)-x_s_2(1)).^2+(X_multi_i(3,:)-x_s_2(2)).^2+(X_multi_i(5,:)-x_s_2(3)).^2);
%         average_dist_sensor2=average_dist_sensor2+dist_sensor;
% 
%     plot(1:Nsteps,dist_sensor,'-r','Linewidth',1.3)
%     
% end
% grid on
% %axis equal
% hold off
% xlabel('Time step')
% ylabel('Distance (m)')
% legend('Sensor 1','Sensor 2')
% Plot average distance to sensors
% average_dist_sensor1=average_dist_sensor1/N_trajectories;
% average_dist_sensor2=average_dist_sensor2/N_trajectories;
% 
% figure(4)
% clf
% plot(1:Nsteps,average_dist_sensor1,'-b','Linewidth',1.3)
% hold on
% plot(1:Nsteps,average_dist_sensor2,'-r','Linewidth',1.3)
% hold off
% grid on
% xlabel('Time step')
% ylabel('Average distance (m)')
% legend('Sensor 1','Sensor 2')
% axis([1 Nsteps 0 600])








