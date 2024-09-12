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



% Dimensionless parameters 
b=cell(n_gene,1);
b{1}=R_constants(1,2)/R_constants(1,3);
     
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


% Initial conditions (Gaussian density function)
PX      = params.P_int;


% Computation of characteristics curves
xbar=cell(n_gene,1);
xbarlim=cell(n_gene,1);
for i=1:n_gene
    xbar{i}=x{i}*exp(deltat*R_constants(i,4));
    xbarlim{i} = find(xbar{i}>=x{i}(end));
end

% Caracteristics grid
Xbargrid=cell(n_gene,1);
[Xbargrid{1:n_gene}] = ndgrid(xbar{1:n_gene});


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
PX_sol      = cell(nt,1);
du1Lix=cell(nt,1);
PX_sol{1}   = PX;

%derivative of the input function on inducers
cx=cell(n_gene,1);
ducx=cell(n_gene,1); 





rho=(b_r/R_constants(1,1))*u(1);

cx{1}=rho+eps;

drho=(b_r/R_constants(1,1));

ducx{1}=drho;

du1Lix{1}=-(R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},ducx{1}.*e_x{1}.*PX,1)-R_constants(1,1)*ducx{1}.*PX);


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

 


expl_den = 1+(R_constants(1,1)*cx{1} - R_constants(1,4))*deltat; 


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
     Lix =R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},e_x{1}.*cx{1}.*PX,1);
       
    % Explicit method    
    PX = (PX_bar+deltat*Lix)./expl_den; 

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
    

% Saving the solution fo the current time step
PX_sol{j} = PX;


%update of the input function
rho=(b_r/R_constants(1,1))*u(j);

cx{1}=rho+eps;

drho=(b_r/R_constants(1,1));

ducx{1}=drho;

%Gradient term
du1Lix{j}=-(R_constants(1,1)/b{1}*e_lx{1}.*cumtrapz(x{1},ducx{1}.*e_x{1}.*PX,1)-R_constants(1,1)*ducx{1}.*PX);
  

end

end 








