%%%%%%%%%%%%%Script for solving the forward problem%%%%%%%%%%%%%%%%%%%%
function [PX_sol]=adjoint_f(params,p_main,p_target,u)


%Parameters extraction
n_gene=params.n_gene;
Prot_mesh=params.protein.Prot_mesh;
Time_mesh=params.time.Time_mesh;
R_constants=params.R_constants;

%Extra parameters 
eps=params.constants.eps;
b_r=params.constants.b_r;


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
tl     = linspace(t0, tmax, nt);
%fprintf('\n The time discretization is: %g \n',tl(2)-tl(1));

% Initial conditions (Gaussian density function)
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
du1Lix=cell(nt_sol,1);

PX_sol{1}   = PX;

u=flip(u);

%%%%%%%%%%filling gradient of xost function t time  0%%%%%%%%%%%%
cx=cell(n_gene,1);

rho=(b_r/R_constants(1,1))*u(1);

cx{1}=rho+eps;




 


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

% Other time independent functions
sumkmcx=0;

for i=1:n_gene
    sumkmcx = sumkmcx + R_constants(i,1)*cx{i};
   
end

expl_den  = 1+(sumkmcx)*deltat;


%%%%%%%%Solving the system for t+dt%%%%%%%%%%%%%

    % PX_bar construction using interpolation 
    if n_gene==1
        PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
        PX_bar(isnan(PX_bar)==1)=0;
    elseif n_gene>1
        PX_bar = interpn(Xgrid{1:n_gene},PX,Xbargrid{1:n_gene});
        PX_bar(isnan(PX_bar)==1)=0;
    end

     
     % Integral term computation by numerical integration
    Lix=0;

    for i=1:n_gene

       Lix = Lix + R_constants(i,1)/b{i}*cx{i}.*e_x{i}.*(-flip(cumtrapz(flip(x{i}),flip(e_lx{i}.*PX,i),i),i));

    end
    
       
    % Explicit method    
    PX = (PX_bar+deltat*(Lix))./expl_den; 
    
    % Zero boundary condition
    CFaux=cell(n_gene,1);
    for i=1:n_gene
        CFaux{i}=':';
    end
    for i=1:n_gene
        CF=CFaux;
        CF{i}=iN{i};
        PX(CF{1:n_gene})=zeros(size(PX(CF{1:n_gene})));
    end

        

       
       %PX=PX./norm(PX);

       %PX=PX./(-(trapz(x{1},PX)));

      
    
    % Saving the solution fo the current time step
    

PX_sol{j} = PX;

rho=(b_r/R_constants(1,1))*u(j);

cx{1}=rho+eps;







   
end


%flip the solution to obtain the backward adjoint 
PX_sol=flip(PX_sol);
end 







