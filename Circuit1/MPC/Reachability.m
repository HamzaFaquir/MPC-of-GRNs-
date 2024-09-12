clear all 
tic
%%%%%%%%% MAIN MPC PROBLEM%%%%%%%%%%%%%%%%%%
n_gene=1;
ti=0;
t_f=100;
tf_max=250;


nt_solve=500; %time disrectization of time windows 
n_x=500;
it_x=5000;
Prot_mesh=[0.000000 n_x it_x; 0.000000 n_x it_x];


deltat = (t_f-t_i)/nt_solve;

iN=cell(n_gene,1);
x=cell(n_gene,1);
for i=1:n_gene
    iN{i} = Prot_mesh(i,3) + 1;
    x{i} = linspace(Prot_mesh(i,1),Prot_mesh(i,2), iN{i});
end


% Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});

%%%%%%%%%%%%%%%%%%%%%%%%%%Target function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_gene = length(Xgrid);


%initial condition 
x0g=0*ones(1,n_gene);  % Mean of the Gaussian density
sigmag=10*ones(1,n_gene); % Standard deviation of the Gaussian density 

gausker=0;
for i=1:n_gene
    gausker = gausker+((Xgrid{i}-x0g(i))/sigmag(i)).^2;
end
PX0un=exp(-gausker/2);

% Norm of the initial condition
auxnor0=PX0un;
for i=1:n_gene
    auxnor0 = trapz(x{i},auxnor0);
end

% Normalization of the initial condition to be a density function
IC=PX0un/auxnor0;
P_k= IC;





mu=100;

PDF = normpdf(x{1}, mu, 30)+ normpdf(x{1}, 25, 15); % take the combination of two normal distrubtions
P = PDF/trapz(x{1},PDF); % normalise it, i.e. sum equals to one.
P=P';
%{
P=load('P_target.mat'); 
P=P.PP;
%}
%%%%GRadient descent optimization%%%%%%%%%%%%%%%%%%%%%%%%%%

tol=10^-2;
n_op=25;
u_= 30;
alpha_step=1000;
u=ones(nt_solve,1)*u_;

[J,P_main, u]=optimal_control_m(ti,t_f,u,P_k,P,Prot_mesh,n_op,nt_solve);

P_optimal=P_main;
u_optimal=u;

test_tol=cost_f(P_main,P,t_f,x);
test_tol_1=test_tol;

J_final=[test_tol_1];
Tf=[t_f];




while t_f<=tf_max && test_tol>tol*test_tol_1

t_f=t_f+10;

nt_solve=round((t_f-ti)/deltat); %time disrectization of time windows 

u=ones(nt_solve,1)*u_;

[J,P_main, u]=optimal_control_m(ti,t_f,u,P_k,P,Prot_mesh,n_op,nt_solve);

test_tol=cost_f(P_main,P,t_f,x)

if sum(J_final>test_tol)==size()
    P_optimal=P_main;
    u_optimal=u;
end 

J_final(end+1)=test_tol;

Tf(end+1)=t_f;


end 
 
fprintf('The optimal value is: %g',min(J_final));
fprintf('\n The minimum reachibility time is: %g ',Tf(J_final==min(J_final)));
 


 
%%%%%%%%%%%%%%Plot results%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=linspace(ti,Tf(J_final==min(J_final)),round((Tf(J_final==min(J_final))-ti)/deltat));




figure

hold on
plot(x{1},P_optimal{end},'-','LineWidth',1.5)
plot(x{1},P,'--')
xlabel('Protein 1')
ylabel('Probability')


hold off

 

%}

%%%%%%%%%%%inducers in time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure

plot(T,u_optimal,'-','LineWidth',1.5)

xlabel('time')
ylabel('u')


if n_op>3
figure

plot(Tf,J_final,'-','LineWidth',1.5)

xlabel('optimisation iterations')
ylabel('Cost function')
end 

toc 

