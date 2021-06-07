%Gaussian tracking with two sensors measuring direction-of-arrival in 3D
%Measurements are Kent distributed. The implemented filter is the iterated posterior
%linearisation filter explained in

% A. F. García-Fernández, S. Maskell, P. Horridge, J. Ralph, “Gaussian tracking with Kent-distributed direction-of-arrival measurements,” IEEE Transactions on Vehicular Technology, 2021.

%This code provides the sigma-point implementation using the unscented
%transform


%Author: Ángel F. García-Fernández

clear
randn('seed',9)
rand('seed',9)
addpath('3D-Simulation-Visualization-master')


Scenario_bearings_3D_Kent;

%Parametros sigma points
Nx=6;
W0=1/3;
Nz=6;% Two 3D measurements for sensor

Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)];
error_square_t_tot=zeros(1,Nsteps);


%Calculation of normalisation constants and mean and covariance when the
%distribution points at [1,0,0]. When we receive the measurement, we just
%need the rotation
log_c=c_Kent(kappa,beta);
log_ck=ck_Kent(kappa,beta);
log_ckk=ckk_Kent(kappa,beta,log_ck);
log_cb=cb_Kent(kappa,beta);
ck_c=exp(log_ck-log_c);  %ck divided by c
ckk_c=exp(log_ckk-log_c);  %ckk divided by c
cb_c=exp(log_cb-log_c);

mean_Kent_100=[ck_c,0,0]';
cov_Kent_100=diag([ckk_c-ck_c^2,(1-ckk_c+cb_c)/2,(1-ckk_c-cb_c)/2]);

if(cov_Kent_100(1)<0)
    cov_Kent_100(1)=0;
end

if(measurement_kent==0)
    R=diag([sigma2_bear,sigma2_el,sigma2_bear,sigma2_el]);
end



nees_t_tot=zeros(1,Nsteps);
error_square_t_series=zeros(Nmc,Nsteps);


Nit_iplf=5; %Number of iterations in the IPLF

mean_ini=x0;

A_l=zeros(6);
b_l=zeros(6,1);
Omega_l=zeros(6);

randn('seed',9)
rand('seed',9)


%Generation of measurements
z_real_t_mc=zeros(Nmc,Nz,Nsteps);
for k=1:Nsteps
    for i=1:N_trajectories
        X_multi=squeeze(X_multi_series(i,:,:));
        state_k=X_multi(:,k);
        %Mode sensor one
        dist_vect1=[state_k(1)-x_s_1(1);state_k(3)-x_s_1(2);state_k(5)-x_s_1(3)];
        mode_z1=dist_vect1/sqrt(sum(dist_vect1.^2));
        
        %Mode sensor two
        dist_vect2=[state_k(1)-x_s_2(1);state_k(3)-x_s_2(2);state_k(5)-x_s_2(3)];
        mode_z2=dist_vect2/sqrt(sum(dist_vect2.^2));
        
        %Measurement generation
        
        if(measurement_kent==1)
           
            %Option 1: Measurements are Kent-distributed
             z_real_1=Random_FB5_Kent(kappa,beta,[1,0,0],pi/2,Nmc_trajectory);
             z_real_1=z_real_1';
             
             z_real_1=RotationKentmeasurements(z_real_1,state_k,x_s_1);  
             
             z_real_2=Random_FB5_Kent(kappa,beta,[1,0,0],pi/2,Nmc_trajectory);
             z_real_2=z_real_2';
             z_real_2=RotationKentmeasurements(z_real_2,state_k,x_s_2);

             z_real=[z_real_1;z_real_2]';
             
            
        else
            %Option 2: Measurements are Gaussian distributed
            z_real_1=[atan2(mode_z1(2),mode_z1(1));asin(mode_z1(3))]+sqrt(R(1:2,1:2))*randn(2,Nmc_trajectory);
            z_real_2=[atan2(mode_z2(2),mode_z2(1));asin(mode_z2(3))]+sqrt(R(1:2,1:2))*randn(2,Nmc_trajectory);
            
            z_real=[cos(z_real_1(1,:)).*cos(z_real_1(2,:));...
                sin(z_real_1(1,:)).*cos(z_real_1(2,:));...
                sin(z_real_1(2,:));...
                cos(z_real_2(1,:)).*cos(z_real_2(2,:));...
                sin(z_real_2(1,:)).*cos(z_real_2(2,:));...
                sin(z_real_2(2,:))]';
        end
        z_real_t_mc((i-1)*Nmc_trajectory+1:i*Nmc_trajectory,:,k)=z_real;
        
        
    end
end


randn('seed',9)
rand('seed',9)

%Number of Monte Carlo runs

for i=1:Nmc
    
    
    
    n_trajectory=fix((i-1)/Nmc_trajectory)+1;
    X_multi=squeeze(X_multi_series(n_trajectory,:,:));
    
    meank=mean_ini;
    Pk=P_ini;
    
    
    error_square_t=zeros(1,Nsteps);
    error_square_t_smoothing=zeros(1,Nsteps);
    
    nees_t=zeros(1,Nsteps);
    
    meank_t=zeros(Nx,Nsteps);
    Pk_t=zeros(Nx,Nx,Nsteps);
    
    
   %Linearisation parameters 
    A_m=zeros(Nz,Nx,Nsteps);
    b_m=zeros(Nz,Nsteps);
    Omega_m=zeros(Nz,Nz,Nsteps);
    

    
    %Measurements
    z_real_t=squeeze(z_real_t_mc(i,:,:));
    
    %VMF with sigma points
    tic
    for k=1:Nsteps
        
        z_real=z_real_t(:,k);
        
        
        %Iterated statistical linear regressions (Nit_iplf steps)
        
        meank_j=meank;
        Pk_j=Pk;
        
        for p=1:Nit_iplf
            
            %Sensor 1
            [A_l1,b_l1,Omega_l1]=SLR_measurement_Kent(meank_j,Pk_j,weights,W0,Nx,mean_Kent_100,cov_Kent_100,x_s_1);
            %Sensor 2
            [A_l2,b_l2,Omega_l2]=SLR_measurement_Kent(meank_j,Pk_j,weights,W0,Nx,mean_Kent_100,cov_Kent_100,x_s_2);
            
            A_l(1:3,:)=A_l1;
            A_l(4:6,:)=A_l2;
            b_l(1:3)=b_l1;
            b_l(4:6)=b_l2;
            Omega_l(1:3,1:3)=Omega_l1;
            Omega_l(4:6,4:6)=Omega_l2;
            
            %KF update (R=0)
            [mean_ukf_act,var_ukf_act]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,0,z_real);
            
            meank_j=mean_ukf_act;
            Pk_j=var_ukf_act;
            
        end
        
        
        A_m(:,:,k)=A_l;
        b_m(:,k)=b_l;
        Omega_m(:,:,k)=Omega_l;
        
        meank_t(:,k)=mean_ukf_act;
        Pk_t(:,:,k)=var_ukf_act;
        
        
        error_square_t(k)=(mean_ukf_act(1)-X_multi(1,k))^2+(mean_ukf_act(3)-X_multi(3,k))^2+(mean_ukf_act(5)-X_multi(5,k))^2;
        pos_error=mean_ukf_act([1,3,5])-X_multi([1,3,5],k);
        var_pos_act=var_ukf_act([1 3 5],[1 3 5]);
        %nees_t(k)=state_error'*inv(var_mp_tukf_act)*state_error;
        nees_t(k)=pos_error'/var_pos_act*pos_error;
       
        
        
        %Prediction
        
        meank=F*mean_ukf_act;
        Pk=F*var_ukf_act*F'+Q;
        Pk=(Pk+Pk')/2;
        
     
        
        
    end
    
    t=toc;

    error_square_t_tot=error_square_t_tot+error_square_t;  
    
    
    nees_t_tot=nees_t_tot+nees_t;
    
    display(['Completed iteration nº ', num2str(i),' time ', num2str(t), ' sec'])
end

error_square_t_tot=error_square_t_tot/Nmc;
error_sqrt_tot=sqrt(sum(error_square_t_tot)/(Nsteps))
nees_t_tot=nees_t_tot/Nmc;





figure(2)
plot(1:Nsteps,sqrt(error_square_t_tot),'b','Linewidth',1.3)
ylabel ('RMS position error (m)')
xlabel('Time step')
grid on
legend('Kent IPLF')


figure(3)
plot(nees_t_tot,'Linewidth',1.3)
grid on
ylabel ('NEES (Kent IPLF)')
xlabel('Time step')
