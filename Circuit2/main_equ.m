%%%%%%%%%%%%%Script for solving the forward problem%%%%%%%%%%%%%%%%%%%%
function [du1Lix, PX_sol]=main_equ(params,u)


%Parameters extraction
n_gene=params.n_gene;
Prot_mesh=params.protein.Prot_mesh;
Time_mesh=params.time.Time_mesh;
R_constants=params.R_constants;

%Extra parameters 
eps=params.constants.eps;
b_r=params.constants.b_r;
k_1=params.constants.k_1;



%%%


% Dimensionless parameters 
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
deltat = (tmax-t0)/nt;


% Initial conditions (Gaussian density function)
PX      = params.P_int;


% Computation of characteristics curves
xbar=cell(n_gene,1);
xbar{1}=x{1}*exp(deltat*params.R_constants(1,4));

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

%Input function and its derivative with respect to inducers
cx=cell(n_gene,1);
ducx=cell(n_gene,1); 

rho=(b_r/k_1)*u(1);
cx{1}=rho+eps;

drho=(b_r/k_1);
ducx{1}=drho;

du1Lix{1}=-(R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},ducx{1}.*e_x{1}.*PX,1)-R_constants(1,1)*ducx{1}.*PX);


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

%%%%%%%%Solving the system for t+dt%%%%%%%%%%%%%
sumkmcx =params.R_constants(1,1)*cx{1};
sumprotdeg =params.R_constants(1,4);
expl_den = 1+(sumkmcx - sumprotdeg)*deltat;

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

% Saving the solution fo the current time step
PX_sol{j} = PX;


%update of the input function
rho=(b_r/k_1)*u(j);

cx{1}=rho+eps;

drho=(b_r/k_1);

ducx{1}=drho;

%Gradient term
du1Lix{j}=-(R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},ducx{1}.*e_x{1}.*PX,1)-R_constants(1,1)*ducx{1}.*PX);
  

end

end 








