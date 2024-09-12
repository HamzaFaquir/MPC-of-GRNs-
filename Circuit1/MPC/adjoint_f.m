%%%%%%%%%%%%%Script for solving the forward problem%%%%%%%%%%%%%%%%%%%%
function [PX_sol]=adjoint_f(params,p_main,p_target,u)

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
    %fprintf('\n The discretization step in the dimension %d is: %g \n',i,x{i}(2)-x{i}(1));
end

% Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});


% Time definition
t0     = Time_mesh(1);
tmax   = Time_mesh(2);
nt     = Time_mesh(3);
deltat = (tmax-t0)/(nt-1);

% Computation of characteristics curves
xbar=cell(n_gene,1);
xbarlim=cell(n_gene,1);
for i=1:n_gene
    xbar{i}=x{i}*exp(-deltat*R_constants(i,4));
    xbarlim{i} = find(xbar{i}>=x{i}(end));
    %fprintf('\nThe length of X%g is: %g and there are %g points of Xbar biggest than Xmax \n',i,length(x{i}),length(xbarlim{i}));
end
% Caracteristics grid
Xbargrid=cell(n_gene,1);
[Xbargrid{1:n_gene}] = ndgrid(xbar{1:n_gene});
 
% Initialization 
PX      =-(p_main{end}-p_target);



% Time Independent functions
e_x=cell(n_gene,1);
e_lx=cell(n_gene,1);
for i=1:n_gene
    e_x{i}=exp(Xgrid{i}/b{i});
    e_lx{i}=exp(-Xgrid{i}/b{i});
    e_x{i}(isinf(e_x{i})==1)=realmax;
    e_lx{i}(isinf(e_lx{i})==1)=realmax;
end



% Saving result for simulated times 
nt_sol = Time_mesh(3);
PX_sol      = cell(nt_sol,1);

PX_sol{1}   = PX;

u=flip(u);

%%%%%%%%%%fCalculation initial input function%%%%%%%%%%%%
cx=cell(n_gene,1);
inducer=(1/(1+(u(1)/K_u)^n_u));

rho=1./(1+((inducer*Xgrid{1})./K_x).^n_x);
rho(isnan(rho)==1)=0;
cx{1}=1-rho+rho*eps;




 


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

% Other time independent functions
expl_den  = 1+(R_constants(1,1)*cx{1})*deltat;


%%%%%%%%Solving the system for t+dt%%%%%%%%%%%%%

% PX_bar construction using interpolation 
PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
PX_bar(isnan(PX_bar)==1)=0;



% Integral term computation by numerical integration
Lix =R_constants(1,1)/b{1}*cx{1}.*e_x{1}.*(-flip(cumtrapz(flip(x{1}),flip(e_lx{1}.*PX,1),1),1));

% Explicit method    
PX = (PX_bar+deltat*(Lix-(p_main{nt-j+1}-p_target)))./expl_den; 

% Zero boundary condition
CFaux=cell(n_gene,1);
CFaux{1}=':';
CF=CFaux;
CF{1}=iN{1};
PX(CF{1:n_gene})=zeros(size(PX(CF{1:n_gene})));
  
PX_sol{j} = PX;


inducer=(1/(1+(u(1)/K_u)^n_u));

rho=1./(1+((inducer*Xgrid{1})./K_x).^n_x);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;
 
end

%flip the solution to obtain the backward adjoint 
PX_sol=flip(PX_sol);
end 







