clear all 
tic
%%%%%%%%% MAIN MPC PROBLEM%%%%%%%%%%%%%%%%%%
n_gene=1;


%Reachability time mesh 
t_int=0;
t_f=100;  %First final time for reachability test 
Reachability_step=50;% step to update the final time 
tf_max=250;  %Maximum final time for reachability test 
nt_solve=1000; 
deltat = (t_f-t_int)/(nt_solve-1);
Time_mesh=[t_int t_f nt_solve];

%Protein mesh
n_x=500;
it_x=5000;
Prot_mesh=[0.000000 n_x it_x; 0.000000 n_x it_x];

%Parameter defintions%%%%%%%%%%%%%%%%%%%%%%%

%Model simulation parameters 
params.time = struct('Time_mesh',Time_mesh);
params.protein = struct('Prot_mesh', Prot_mesh);
params.constants=struct('K_x',80,'n_x',4,'n_u',3,'K_u',200,'eps',0.1);
params.R_constants = [10.0000, 1.2500, 0.2500, 0.1750];
params.n_gene=n_gene;


%Optimization parameters 
params.opt.alpha0=10^12;
params.opt.alpha_kmax=60;
params.opt.number_operations=25;
params.opt.NCG_tol=10^(-5);
params.opt.Inducer_upbo=300;
params.opt.Inducer_lobo=0;
params.opt.step_search=struct('gamma',10^(-7),'beta',0.7);



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
P_0= IC;
params.P_int=P_0;


%%%colours for plot%%%%
str1 = '#77AC30';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;

str2 = '#A2142F';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;

str3 = '#EDB120';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;



%%%Setting the target to test reachability 
mu=130;

PDF = normpdf(x{1}, mu, 15)+ normpdf(x{1}, 25, 15); % take the combination of two normal distrubtions
P = PDF/trapz(x{1},PDF); % normalise.
P=P';
params.opt.P_target=P;

%%%%GRadient descent optimization%%%%%%%%%%%%%%%%%%%%%%%%%%


tol=10^-2;
n_op=25;
u_= 30;
u=ones(nt_solve,1)*u_;

%%%First simulation of the first final time
[J,P_main, u]=optimal_control_m(params,u);

P_optimal=P_main;
u_optimal=u;

P_opt=P_main{end};

%{
filename = sprintf('P_opt_tf%d.mat',t_f);
save(filename,'P_opt')

filename = sprintf('u_opt_tf%d.mat',t_f);
save(filename,'u_optimal')
%}

test_tol=cost_f(P_main,P,t_f,x);
test_tol_1=test_tol;
test_optimal=test_tol_1;

plot(x{1},P_main{end},'-','Color',color1,'LineWidth',2.5)
hold on 
plot(x{1},P,'--','LineWidth',2.5)
xlabel('Protein 1')
ylabel('Probability')
 



%functions to store the cost function for its final time with the
%corresponding final time

J_final=[test_tol_1];
Tf=[t_f];



filename1 = sprintf('reach_Pd_mu=%d.eps',mu);
filename2 = sprintf('reach_u_mu=%d.eps',mu);

filename3 = sprintf('reach_J_mu=%d.eps',mu);

while t_f<=tf_max 
   
fprintf('\n The current time is: %g',t_f);
fprintf('\n The current optimal value is: %g',min(J_final));
fprintf('\n The current optimal time is: %g \n ',Tf(J_final==min(J_final)));

if test_tol>tol*test_tol_1

t_f=t_f+Reachability_step;

nt_solve=round(1+(t_f-t_int)/deltat); %time disrectization of time windows 
params.time = struct('Time_mesh',[t_int t_f nt_solve]);

u=ones(nt_solve,1)*u_;

[J,P_main, u]=optimal_control_m(params,u);

test_tol=cost_f(P_main,P,t_f,x);

P_opt=P_main{end};
%{
filename = sprintf('P_opt_tf%d.mat',t_f);
save(filename,'P_opt')

filename = sprintf('u_opt_tf%d.mat',t_f);
save(filename,'u_optimal')
%}

if test_optimal>test_tol
    test_optimal=test_tol;
    P_optimal=P_main;
    u_optimal=u;
end 

J_final(end+1)=test_tol;

Tf(end+1)=t_f;

else 
    disp('Tolerance reached')
end 

end 
 
fprintf('The optimal value is: %g',min(J_final));
fprintf('\n The minimum reachibility time is: %g ',Tf(J_final==min(J_final)));
 


 
%%%%%%%%%%%%%%Plot results%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=linspace(Time_mesh(1),Tf(J_final==min(J_final)),round((Tf(J_final==min(J_final))-Time_mesh(1))/deltat)+1);

P_reachable=P_optimal{end};
filename = sprintf('P_reachable_mu%d.mat',mu);
save(filename,'P_reachable')

filename = sprintf('u_reachable_mu%d.mat',mu);
save(filename,'u_optimal')


filename = sprintf('J_reachable_mu%d.mat',mu);
save(filename,'J_final')

%%%%plot results%%%%%%%%%%%%
str1 = '#77AC30';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;

str2 = '#A2142F';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;

str3 = '#EDB120';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;


figure1=figure;

hold on
plot(x{1},P_optimal{end},'-','Color',color1,'LineWidth',2.5)
plot(x{1},P,'--','LineWidth',2.5)
xlabel('Protein 1')
ylabel('Probability')

legend('P_reach','Target')

%saveas(figure1,filename1)  % here you save the figure

hold off

 

figure2=figure;
hold on
plot(T,u_optimal,'-','Color',color2,'LineWidth',2.5)

xlabel('time')
ylabel('u')

%saveas(figure2,filename2)  % here you save the figure
hold off 

if n_op>3
figure3=figure;
hold on 
plot(Tf,J_final,'-','Color',color3,'LineWidth',2.5)

xlabel('Final times')
ylabel('Cost function')

hold off 
%saveas(figure3,filename3)  % here you save the figure

end 

toc 

