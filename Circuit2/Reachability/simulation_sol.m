clear all 

%%%%%%%%% MAIN MPC PROBLEM%%%%%%%%%%%%%%%%%%
n_gene=1;
nt_solve=1000; %time disrectization of time windows 
t_int=0;
t_f=100;



n_x=100;
it_x=1000;
Prot_mesh=[0.000000 n_x it_x];

Time_mesh=[t_int t_f nt_solve];

T=linspace(t_int,t_f, nt_solve);



iN=cell(n_gene,1);
x=cell(n_gene,1);
for i=1:n_gene
    iN{i} = Prot_mesh(i,3) + 1;
    x{i} = linspace(Prot_mesh(i,1),Prot_mesh(i,2), iN{i});
end

%%%
%Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});

%%%%%%%%%%%%%%%%%%%%%%%%%%Target function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_gene = length(Xgrid);


iN=cell(n_gene,1);
x=cell(n_gene,1);
for i=1:n_gene
    iN{i} = Prot_mesh(i,3) + 1;
    x{i} = linspace(Prot_mesh(i,1),Prot_mesh(i,2), iN{i});
end

%%%
%Protein "spatial" Mesh
Xgrid=cell(n_gene,1);
[Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});

%%%%%%%%%%%%%%%%%%%%%%%%%%Target function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initial condition 
x0g=0*ones(1,n_gene);  % Mean of the Gaussian density
sigmag=5*ones(1,n_gene); % Standard deviation of the Gaussian density 

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


%Parameter defintions
params.time = struct('Time_mesh',Time_mesh);
params.protein = struct('Prot_mesh', Prot_mesh);
params.constants=struct('eps',0.5,'b_r',0.0965,'n_gene',1);
params.R_constants = [0.0048 ,0.0116,0.0048,0.0016];
params.P_int=P_0;
params.n_gene=n_gene;


%PDF = 0.1*normpdf(x{1}, 50, 10);% + 30*normpdf(x{1}, 350, 25); % take the combination of two normal distrubtions
%P = PDF;%/sum(PDF); % normalise it, i.e. sum equals to one.
%P=P';

u=ones(1,nt_solve)*20;
[du1Lix,P_main]=main_equ(params,u); %solution of the main problem 

%trapz(x{1},P_main{end})
figure

hold on
plot(x{1},P_main{end},'-','LineWidth',1.5)
%plot(x{1},P,'--','LineWidth',1.5)
xlabel('Protein 1')
ylabel('Probability')

hold off



 
