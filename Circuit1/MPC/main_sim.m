%%%%%%%%%%%%%Script for solving the forward problem%%%%%%%%%%%%%%%%%%%%
function [PX]=main_sim(t_int,t_f,n_tt,Prot_mesh,u,PX0un)

%production rates mrna
a=20;

d_r=0.5;
eps=0.1;

%protein
b_p=25*0.1;
d_p=0.35;

H=4;
K=80;










n_gene=1;
%Prot_mesh=[0.000000 n_x it_x; 0.000000 n_x it_x];
Time_mesh=[t_int t_f n_tt 1];
R_constants=[a b_p d_r d_p]./2;



% Numerical solution of the Friedman equation in general dimension with 
% the semilagrangian method
%%%


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
nt     = Time_mesh(3)*Time_mesh(4);
deltat = (tmax-t0)/nt;
tl     = linspace(t0, tmax, nt);
%fprintf('\n The time discretization is: %g \n',tl(2)-tl(1));

% Initial conditions (Gaussian density function)

% Normalized initial condition
auxnor0=PX0un;
for i=1:n_gene
    auxnor0 = trapz(x{i},auxnor0);
end
% initial condition normalized
PX0=PX0un/auxnor0;

PX      = PX0;


% Computation of characteristics curves
xbar=cell(n_gene,1);
xbarlim=cell(n_gene,1);
for i=1:n_gene
    xbar{i}=x{i}*exp(deltat*R_constants(i,4));
    xbarlim{i} = find(xbar{i}>=x{i}(end));
    %fprintf('\nThe length of X%g is: %g and there are %g points of Xbar biggest than Xmax \n',i,length(x{i}),length(xbarlim{i}));
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







%%%%%%%%%%filling gradient of xost function t time  0%%%%%%%%%%%%
cx=cell(n_gene,1);

c=3;
L=200;

inducer=1/(1+(u(1)/L)^c);

%rho=(inducer*Xgrid{1}.^H)./(inducer*Xgrid{1}.^H+K.^H);
rho=1./(1+((inducer*Xgrid{1})./K).^H);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;


%%%%%%%%%%%%%Time iteration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:nt

    % Other time independent functions
sumkmcx=0;

for i=1:n_gene
    sumkmcx = sumkmcx + R_constants(i,1)*cx{i};%.*(1-exp(Xgrid{i}-20*Prot_mesh(i,2)/b{i}));
   
end
sumprotdeg=0;
for i=1:n_gene
    sumprotdeg = sumprotdeg + R_constants(i,4);
end

expl_den = 1+(sumkmcx - sumprotdeg)*deltat; 


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
        Lix = Lix + R_constants(i,1)/b{i}*e_lx{i}.*cumtrapz(x{i},e_x{i}.*cx{i}.*PX,i);
        
    end
    %saving the gradient term  for the optimization
    
       
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

inducer=1/(1+(u(j)/L)^c);

rho=1./(1+((inducer*Xgrid{1})./K).^H);

rho(isnan(rho)==1)=0;

cx{1}=1-rho+rho*eps;



end

end 








