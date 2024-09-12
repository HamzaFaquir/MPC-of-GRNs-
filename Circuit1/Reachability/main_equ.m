%%%%%%%%%%%%%Script for solving the forward problem%%%%%%%%%%%%%%%%%%%%
function [du1Lix, PX_sol]=main_equ(params,u)

%Extract parameters
n_gene=params.n_gene;
Prot_mesh=params.protein.Prot_mesh;
Time_mesh=params.time.Time_mesh;
R_constants=params.R_constants;

eps=params.constants.eps;
K_x=params.constants.K_x;
n_x=params.constants.n_x;
K_u=params.constants.K_u;
n_u=params.constants.n_u;


% Numerical solution of the Friedman equation in general dimension with 
% the semilagrangian method
%%%


% Dimensionless parameters (gamma1 = 0.01 and gamma2 = 4e-4)
b=cell(n_gene,1);

b{1}=R_constants(1,2)/R_constants(1,3);
     

% Spatial discretization
iN=cell(n_gene,1);
x=cell(n_gene,1);

iN{1} = Prot_mesh(1,3) + 1;
x{1} = linspace(Prot_mesh(1,1),Prot_mesh(1,2), iN{1});
    
% Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});


% Time definition
t0     = Time_mesh(1);
tmax   = Time_mesh(2);
nt     = Time_mesh(3);
deltat = (tmax-t0)/(nt-1);


% Initial conditions
PX= params.P_int;


% Computation of characteristics curves
xbar=cell(n_gene,1);
xbarlim=cell(n_gene,1);
xbar{1}=x{1}*exp(deltat*R_constants(1,4));
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




% Saving result for simulated times 
PX_sol      = cell(nt,1);

du1Lix=cell(nt,1);
PX_sol{1}   = PX;


%%%%%%%%%%Calculation of gradient of inputfunction dudc%%%%%%%%%%%%
cx=cell(n_gene,1);
ducx=cell(n_gene,1); 

inducer=1/(1+(u(1)/K_u)^n_u);

dinducer=((-n_u/K_u)*(u(1)/K_u)^(n_u-1))/((1+(u(1)/K_u)^n_u)^2);

rho=1./(1+((inducer*Xgrid{1})./K_x).^n_x);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;

drho=(-n_x*(dinducer.*Xgrid{1}./K_x).*((inducer*Xgrid{1})./K_x).^(n_x-1))./((1+((inducer*Xgrid{1})./K_x).^n_x).^2);

ducx{1}=drho*(eps-1);

ducx{1}(isnan(ducx{1})==1)=0;

du1Lix{1}=-(R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},ducx{1}.*e_x{1}.*PX,1)-R_constants(1,1)*ducx{1}.*PX);


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

    
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
    

% Saving the solution fo the current time step
PX_sol{j} = PX;

%%%%%%%%%%Calculation of gradient of inputfunction dudc%%%%%%%%%%%%
inducer=1/(1+(u(j)/K_u)^n_u);
dinducer=(-n_u/K_u)*(u(j)/K_u)^(n_u-1)/((1+(u(j)/K_u)^n_u)^2);

rho=1./(1+((inducer*Xgrid{1})./K_x).^n_x);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;

drho=(-n_x*(dinducer.*Xgrid{1}./K_x).*((inducer*Xgrid{1})./K_x).^(n_x-1))./((1+((inducer*Xgrid{1})./K_x).^n_x).^2);

ducx{1}=drho*(eps-1);

ducx{1}(isnan(ducx{1})==1)=0;

du1Lix{j}=-(R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},ducx{1}.*e_x{1}.*PX,1)-R_constants(1,1)*ducx{1}.*PX);

end
end 








