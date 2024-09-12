%%%%%%%%%%%%% Model predictive control of circuit 2%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Code for the in silico MPC of Ciurcuit 2 %%%%%%%%%%%%%%%%%
%%%%%%%%Kinetics of the network: c(x,u)= eps+(b_r/k_19*u) %%%%%%%%%%%%%

close all
clear all

n_gene=1;



t_int=0;
T_sample=20; %sampling time
m=200; %discretization of sample interval
N_max=200; %N_max*T_sample the number of sample times to run mpc
N_ch=1; %control horizon T*N_ch is the time window used for prediction and control 


n_x=500; %protein number
it_x=5000; %discretizatiom number of protein mesh 

Prot_mesh=[0.000000 n_x it_x];
Time_mesh=[t_int T_sample m 1];
deltat = (Time_mesh(2)-Time_mesh(1))/(Time_mesh(3)-1);


%%%%%Parameter definition%%%%%%%%
%Model simulation parameters 
params.time = struct('Time_mesh',Time_mesh);
params.protein = struct('Prot_mesh', Prot_mesh);
params.constants=struct('eps',0.5,'b_r',0.0965,'n_gene',1,'k_1',0.0956);
params.R_constants = [0.0048 ,0.0116,0.0048,0.0016];
params.n_gene=n_gene;

%optimal control parameters 
params.opt.alpha0=10^7;
params.opt.alpha_kmax=20;
params.opt.number_operations=25;
params.opt.NCG_tol=10^(-5);
params.opt.Inducer_upbo=250;
params.opt.Inducer_lobo=0;
params.opt.step_search=struct('gamma',(10^(-4)),'beta',0.7);

%Mpc parameters 
params.mpc.N_ch=N_ch;



% Dimensionless parameters 
b=cell(n_gene,1);

b{1}=params.R_constants(1,2)/params.R_constants(1,3);
    

% Spatial discretization
iN=cell(n_gene,1);
x=cell(n_gene,1);

iN{1} = Prot_mesh(1,3) + 1;
x{1} = linspace(Prot_mesh(1,1),Prot_mesh(1,2), iN{1});


% Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});





% Initial conditions (Gaussian density function)
x0g=0*ones(1,n_gene);  % Mean of the Gaussian density
sigmag=10*ones(1,n_gene); % Standard deviation of the Gaussian density 
gausker =((Xgrid{1}-x0g(1))/sigmag(1)).^2;
PX0un=exp(-gausker/2);
%Normalization of the initial condition to be a density function
PX0un=PX0un/trapz(x{1},exp(-gausker/2));
PX=PX0un;
%loading initial conditions
params.P_int=PX;

% Computation of characteristics curves
xbar=cell(n_gene,1);
xbarlim=cell(n_gene,1);
xbar{1}=x{1}*exp(deltat*params.R_constants(1,4));
xbarlim{1} = find(xbar{1}>=x{1}(end));

% Caracteristics grid
Xbargrid=cell(n_gene,1);
[Xbargrid{1:n_gene}] = ndgrid(xbar{1:n_gene});

% Time Independent functions
e_x=cell(n_gene,1);
e_lx=cell(n_gene,1);
e_x{1}=exp(Xgrid{1}/b{1});
e_lx{1}=exp(-Xgrid{1}/b{1});
e_x{1}(isinf(e_x{1})==1)=realmax;
e_lx{1}(isinf(e_lx{1})==1)=realmax;



% Cell to save solutions at each time step 
PX_sol      = cell(m*N_max,1);
PX_sol{1}   = PX;
P_k=PX;


%%%%%%%%%%filling gradiem*N_max of xost function t time  0%%%%%%%%%%%%
cx=cell(n_gene,1);
ducx=cell(n_gene,1); %derivative of the input function on inducers


option=2; %Option 1 to load a reachable target from the folder 2 to define a new target

if option==1
mu=100;
filename = sprintf('P_reachable_mu%d.mat',mu);
P_target=load(filename);
P=P_target.P_reachable;
P=P/trapz(x{1},P);
params.opt.P_target=P;
else 
mu=120;
PDF = 1*normpdf(x{1}, mu, 18);%
P = PDF/trapz(x{1},PDF); % normalise it, i.e. sum equals to one.
P=P';
params.opt.P_target=P;
end 

%Initsialization of the input 
u_sol=zeros(m*N_max,1);

u_1=10;
u_sol(1)=u_1;


cost_function=zeros(m*N_max,1);
cost_function(1)=max(abs(PX_sol{1}-P));




cx{1}=params.constants.eps+(params.constants.b_r/params.constants.k_1)*u_sol(1);
ducx{1}=(params.constants.b_r/params.constants.k_1);


% Create file name variable to animation 
filename = sprintf('Evolution of Pd for T=%d CH=%d mu=%d.gif', T_sample,N_ch,mu);
filename2 = sprintf('Evolution of u for T=%d CH=%d mu=%d.gif', T_sample,N_ch,mu);
filename3 = sprintf('Evolution of J for T=%d CH=%d mu=%d.gif', T_sample,N_ch,mu);

figure

%%%%%%%%%%%%%Running MPC and Simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
for j = 1:m*N_max
drawnow

%%%%Running controller at every sampling Time%%%%%%%%%% 
if mod(j,m)==0 || j==1
        
if j~= m*N_max
        
    params.time = struct('Time_mesh',[T_sample*(j-1) T_sample*(j-1)+T_sample*(N_ch) N_ch*m m]);
    params.P_int=P_k;
    [J,P_main, u_]=optimal_control_mpc(params,ones(N_ch*m,1)*u_1);
    
end 
end  
 
%%%%Control input update%%%%%%%%%%%% 
u_sol(j)=u_(1);


%%Update of Input function%%%%%%%%%%%%%%%%
cx{1}=params.constants.eps+(params.constants.b_r/params.constants.k_1)*u_sol(j);

%%%%%%%%Solving the system for t+dt%%%%%%%%%%%%%
expl_den = 1+(params.R_constants(1,1)*cx{1} - params.R_constants(1,4))*deltat;

% PX_bar construction using im*N_maxerpolation 
PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
PX_bar(isnan(PX_bar)==1)=0;
    
     
Lix =params.R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},e_x{1}.*cx{1}.*PX,1);
    
% Explicit method    
PX = (PX_bar+deltat*Lix)./expl_den; 
    
% Zero boundary condition
CFaux=cell(n_gene,1);
CFaux{1}=':';
CF=CFaux;
CF{1}=iN{1};
PX(CF{1:n_gene})=zeros(size(PX(CF{1:n_gene})));
  
PX=PX/trapz(x{1},PX);    

% Saving the solution fo the current*N_max time step
PX_sol{j} = PX;
P_k=PX;
cost_function(j)=max(abs(PX_sol{j}-P));

%%%%Plot solution current time and animation    
 if mod(j,m)==0 || j==1
     
     plot(linspace(t_int,T_sample*j/m,j),u_sol(1:j,1),'-','LineWidth',1.5)
     xlabel('time')
     ylabel('input')
     xlim([0 T_sample*N_max])
     ylim([0 500])

     title(sprintf('Time: %0.2f a.u, mu: %f', j*deltat,mu),...
    'Interpreter','Latex');
    
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j==1
  
        imwrite(imind,cm,filename2,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
    imwrite(imind,cm,filename2,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end 
    hold off

 
     plot(linspace(t_int,T_sample*j/m,j),cost_function(1:j,1),'-','LineWidth',1.5)
     xlabel('time')
     ylabel('cost function')
     xlim([0 T_sample*N_max])
     ylim([0 max(cost_function)*1.2])

     title(sprintf('Time: %0.2f a.u, mu: %f', j*deltat,mu),...
    'Interpreter','Latex');
    
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j==1
  
        imwrite(imind,cm,filename3,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
    imwrite(imind,cm,filename3,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end 
    hold off

     plot(x{1},PX,'-','LineWidth',1.5)
     hold on 
     plot(x{1},P,'--','LineWidth',1.5)
     xlabel('Protein 1')
     ylabel('Probability')
     title(sprintf('Time: %0.2f a.u, mu: %f', j*deltat,mu),...
    'Interpreter','Latex');

    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j==1
  
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',0.1);
    else
    imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.1);
    end 
    hold off
 end 
  
end


%Saving solutions 
 filename1 = sprintf('P_mpc_T%d_ch%d_mu%d.mat',T_sample,N_ch,mu);
 filename2 = sprintf('u_mpc_T%d_ch%d_mu%d.mat',T_sample,N_ch,mu);
 filename3 = sprintf('cost_mpc_T%d_ch%d_mu%d.mat',T_sample,N_ch,mu);


 save(filename1,'PX_sol');
 save(filename2,'u_sol');
 save(filename3,'cost_function');

%Plot final solution 
figure
hold on
plot(x{1},PX_sol{end},'-','LineWidth',1.5)
plot(x{1},P,'--','LineWidth',1.5)
xlabel('Protein 1')
ylabel('Probability')
hold off

figure
plot(linspace(ti,T_sample*N_max,m*N_max),u_sol,'-','LineWidth',1.5)
xlabel('time')
ylabel('u')


if n_op>3
figure
plot(linspace(ti,T_sample*N_max,m*N_max),cost_function,'-','LineWidth',1.5)
xlabel('optimisation iterations')
ylabel('Cost function')
end 

toc
