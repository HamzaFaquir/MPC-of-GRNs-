clear all 
tic
%%%%%%%%% MAIN MPC PROBLEM%%%%%%%%%%%%%%%%%%
n_gene=1;
ti=0;
tf=200;



nt_solve=2000; %time disrectization of time windows 
n_x=500;
it_x=5000;
Prot_mesh=[0.000000 n_x it_x; 0.000000 n_x it_x];




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



mu=130;

filename = sprintf('P_reachable_mu%d.mat',mu);
P_target=load(filename);
P=P_target.P_reachable;


%{ 
PDF = normpdf(x{1}, mu, 15)+ normpdf(x{1}, 25, 15); % take the combination of two normal distrubtions
P = PDF/trapz(x{1},PDF); % normalise.
P=P';
%}





%%%%GRadient descent optimization%%%%%%%%%%%%%%%%%%%%%%%%%%
n_op=25;
u_= 30;

u=ones(nt_solve,1)*u_;


problem=1;


%Solve optimization problem 
if problem==1

%u=randi(300,nt_solve,1);
[J,P_main, u]=optimal_control_m(ti,tf,u,P_k,P,Prot_mesh,n_op,nt_solve);

elseif problem == 2
 ti=0;
 T_sample=10;

 n=20;


 nt_solve=100; %time disrectization of time windows 
 tf=T_sample*n;

 u=ones(nt_solve*n,1)*u_;

[J,P_main, u]=optimal_control_mpc(ti,tf,u,P_k,P,Prot_mesh,n_op,n,nt_solve);
nt_solve=nt_solve*n; 

else 

[J,P_main, u_]=optimal_control_m_c(ti,tf,u_,P_k,P,Prot_mesh,n_op,nt_solve);

u=ones(nt_solve,1)*u_;

end 
 




 
%%%%%%%%%%%%%%Plot results%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=(find(J~=0, 1, 'last' ));
if k==0
    k=n_op;
end 
T=linspace(ti,tf, nt_solve);
It=linspace(1,n_op, n_op);


str1 = '#77AC30';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;

str2 = '#A2142F';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;

str3 = '#EDB120';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;



figure

hold on
plot(x{1},P_main{end},'-','Color',color1,'LineWidth',2.5)
plot(x{1},P,'--')
legend('optimal P','Target P')
xlabel('Protein 1')
ylabel('Probability')


hold off


 

%}

%%%%%%%%%%%inducers in time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure

plot(T,u,'-','Color',color2,'LineWidth',2.5)

xlabel('time')
ylabel('u')


if n_op>3
figure

plot(It(1:k),J(1:k),'-','Color',color3,'LineWidth',2.5)

xlabel('optimisation iterations')
ylabel('Cost function')
end 

toc 