%%%%%%%%%%%%%Script for solving the backward adjoint problem%%%%%%%%%%%%%%%%%%%%
function [PX_sol]=adjoint_f(params,p_main,p_target,u)


%Parameters extraction
n_gene=params.n_gene;
Prot_mesh=params.protein.Prot_mesh;
Time_mesh=params.time.Time_mesh;
R_constants=params.R_constants;

%Extra parameters 
eps=params.constants.eps;
b_r=params.constants.b_r;
k_1=params.constants.k_1;


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

% Computation of characteristics curves
xbar=cell(n_gene,1);
xbar{1}=x{1}*exp(-deltat*R_constants(1,4));

% Caracteristics grid
Xbargrid=cell(n_gene,1);
[Xbargrid{1:n_gene}] = ndgrid(xbar{1:n_gene});
 
% Initialization 
PX =-(p_main{end}-p_target);



% Time Independent functions
e_x=cell(n_gene,1);
e_lx=cell(n_gene,1);

e_x{1}=exp(Xgrid{1}/b{1});
e_lx{1}=exp(-Xgrid{1}/b{1});
e_x{1}(isinf(e_x{1})==1)=realmax;
e_lx{1}(isinf(e_lx{1})==1)=realmax;




% Saving result for simulated times 
nt_sol = Time_mesh(3);
PX_sol      = cell(nt_sol,1);

PX_sol{1}   = PX;

u=flip(u);

%%%%%%%%%% Input function c(x,u)%%%%%%%%%%%%
cx=cell(n_gene,1);
rho=(b_r/k_1)*u(1);
cx{1}=rho+eps;


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

%%%%%%%%Solving the system for t+dt%%%%%%%%%%%%%    
sumkmcx =  R_constants(1,1)*cx{1};
   
expl_den  = 1+(sumkmcx)*deltat;

PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
PX_bar(isnan(PX_bar)==1)=0;

Lix =  R_constants(1,1)/b{1}*cx{1}.*e_x{1}.*(-flip(cumtrapz(flip(x{1}),flip(e_lx{1}.*PX,1),1),1));    
       
% Explicit method    
PX = (PX_bar+deltat*(Lix))./expl_den; 
    
% Zero boundary condition
CFaux=cell(n_gene,1);
CFaux{1}=':';
CF=CFaux;
CF{1}=iN{1};
PX(CF{1:n_gene})=zeros(size(PX(CF{1:n_gene})));

 
% Saving the solution fo the current time step
PX_sol{j} = PX;

rho=(b_r/k_1)*u(j);
cx{1}=rho+eps;
   
end


%flip the solution to obtain the backward adjoint 
PX_sol=flip(PX_sol);
end 







