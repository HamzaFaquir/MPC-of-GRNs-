%%%%%%%%%%%%%Script for solving the forward problem%%%%%%%%%%%%%%%%%%%%
close all
clear all

n_gene=1;

%Defining main Time and protein mesh paramters 
ti=0;
T_sample=1; %sample time
m=10; %discretization of sample interval
N_max=30; %N_max defined as N_max*m the numver of iteration of mpc
N_ch=4; %control horizon

n_x=500; %protein number
it_x=5000; %discretizatiom number of protein mesh 

Prot_mesh=[0.000000 n_x it_x; 0.000000 n_x it_x];
Time_mesh=[ti T_sample m 1];


%Loading parameters
%Model simulation parameters 
params.time = struct('Time_mesh',Time_mesh);
params.protein = struct('Prot_mesh', Prot_mesh);
params.constants=struct('K_x',80,'n_x',4,'n_u',3,'K_u',200,'eps',0.1);
params.R_constants = [10.0000, 1.2500, 0.2500, 0.1750];
params.n_gene=n_gene;

R_constants=params.R_constants;
eps=params.constants.eps;
K_x=params.constants.K_x;
n_x=params.constants.n_x;
K_u=params.constants.K_u;
n_u=params.constants.n_u;

%optimal control parameters 
params.opt.alpha0=10^12;
params.opt.alpha_kmax=20;
params.opt.number_operations=25;
params.opt.NCG_tol=10^(-8);
params.opt.Inducer_upbo=300;
params.opt.Inducer_lobo=0;
params.opt.step_search=struct('gamma',(10^(-7)),'beta',0.7);

%Mpc parameters 
params.mpc.N_ch=N_ch;


% Dimensionless parameters (gamma1 = 0.01 and gamma2 = 4e-4)
b=cell(n_gene,1);
for i=1:n_gene
    b{i}=R_constants(i,2)/R_constants(i,3);
end      

% Spatial discretization
iN=cell(n_gene,1);
x=cell(n_gene,1);
for i=1:n_gene
    iN{i} = Prot_mesh(i,3) + 1;
    x{i} = linspace(Prot_mesh(i,1),Prot_mesh(i,2), iN{i});
    %fprim*N_maxf('\n The discretization step in the dimension %d is: %g \n',i,x{i}(2)-x{i}(1));
end

% Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});


% Time definition
t0     = Time_mesh(1);
tmax   = Time_mesh(2);
deltat = (tmax-t0)/(Time_mesh(3)-1);




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
A=zeros(N_max,iN{1});

PX_sol{1}   = PX;
A(1,:)=PX;
P_k=PX;

%%%%%%%%%%filling gradiem*N_max of xost function t time  0%%%%%%%%%%%%
cx=cell(n_gene,1);
ducx=cell(n_gene,1); %derivative of the input function on inducers



option=1; %Option 1 to load a reachable target from the folder 2 to define a new target

if option==1
%Load the bimadol distrbution
mu=200; %second mode
filename = sprintf('P_reachable_mu%d.mat',mu);
P_target=load(filename);
P=P_target.P_reachable;
P=P/trapz(x{1},P);
params.opt.P_target=P;
else 
%Create a normal distribution 
mu=120;
PDF = 1*normpdf(x{1}, mu, 18);%
P = PDF/trapz(x{1},PDF); % normalise it, i.e. sum equals to one.
P=P';
params.opt.P_target=P;
end 


u_sol=zeros(m*N_max,1);

u_1=10;
u_sol(1)=u_1;


cost_function=zeros(m*N_max,1);
cost_function(1)=max(abs(PX_sol{1}-P));



inducer=1/(1+(u_sol(1)/K_u)^n_u);

rho=1./(1+((inducer*Xgrid{1})./K_x).^n_x);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;


% Create file name variable
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

u_sol(j)=u_(1);

inducer=1/(1+(u_sol(j)/K_u)^n_u);

rho=1./(1+((inducer*Xgrid{1})./K_x).^n_x);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;



% Other time independent functions
expl_den = 1+(R_constants(1,1)*cx{1} - R_constants(1,4))*deltat; 
PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
PX_bar(isnan(PX_bar)==1)=0;

%%%%%%%%Solving the system for t+dt%%%%%%%%%%%%%

% Integral term computation by numerical integration  
Lix =R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},e_x{1}.*cx{1}.*PX,1);
            
    
% Explicit method    
PX = (PX_bar+deltat*Lix)./expl_den; 

% Zero boundary condition
CFaux=cell(n_gene,1);
CFaux{1}=':';
CF=CFaux;
CF{1}=iN{1};
PX(CF{1:n_gene})=zeros(size(PX(CF{1:n_gene})));
        
PX=PX/trapz(x{1},PX);

% Saving the solution fo the currem*N_max time step
PX_sol{j} = PX;
P_k=PX;
cost_function(j)=max(abs(PX_sol{j}-P));

  

%Plot current solution and generate animation of MPC
 if mod(j,m)==0 || j==1
     
     plot(linspace(ti,T_sample*j/m,j),u_sol(1:j,1),'-','LineWidth',1.5)
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
        'DelayTime',0.5);
    else
    imwrite(imind,cm,filename2,'gif','WriteMode','append',...
        'DelayTime',0.5);
    end 
    hold off

   
    
     plot(linspace(ti,T_sample*j/m,j),cost_function(1:j,1),'-','LineWidth',1.5)
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
        'DelayTime',0.5);
    else
    imwrite(imind,cm,filename3,'gif','WriteMode','append',...
        'DelayTime',0.5);
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
        'DelayTime',0.5);
    else
    imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.5);
    end 
    hold off


 end  


end

 filename1 = sprintf('P_mpc_T%d_ch%d_mu%d.mat',T_sample,N_ch,mu);
 filename2 = sprintf('u_mpc_T%d_ch%d_mu%d.mat',T_sample,N_ch,mu);
 filename3 = sprintf('cost_mpc_T%d_ch%d_mu%d.mat',T_sample,N_ch,mu);


 save(filename1,'PX_sol');
 save(filename2,'u_sol');
 save(filename3,'cost_function');


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


if params.opt.number_operations>3
figure

plot(linspace(ti,T_sample*N_max,m*N_max),cost_function,'-','LineWidth',1.5)

xlabel('optimisation iterations')
ylabel('Cost function')
end 
toc 

 str1 = '#77AC30';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;

str2 = '#A2142F';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;

str3 = '#EDB120';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;


str4 = '#0072BD';
color4 = sscanf(str4(2:end),'%2x%2x%2x',[1 3])/255;

figure 
[X,Y] = meshgrid(linspace(ti,T_sample*N_max,m*N_max),x{1});
plot3(X,Y,PX_sol{end},'Color',color3,LineWidth=2.5)
hold on
%plot3(X,Y,[P_int P_int P2 P2 P3 P3],'--','Color',color1,LineWidth=2.5)
xlabel('Time')
ylabel('Protein')
zlabel('Probability')
grid on 
view(60,60);




