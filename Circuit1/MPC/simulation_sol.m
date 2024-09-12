clear all 

%%%%%%%%% MAIN MPC PROBLEM%%%%%%%%%%%%%%%%%%
n_gene=1;
nt_solve=3000; %time disrectization of time windows 
ti=0;
tf=300;
m=10;


n_x=500;
it_x=5000;
Prot_mesh=[0.000000 n_x it_x];

Time_mesh=[ti tf nt_solve 1];

T=linspace(ti,tf, nt_solve);



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
P_k= IC;


%Model simulation parameters 
params.time = struct('Time_mesh',Time_mesh);
params.protein = struct('Prot_mesh', Prot_mesh);
params.constants=struct('K_x',80,'n_x',4,'n_u',3,'K_u',200,'eps',0.1);
params.R_constants = [10.0000, 1.2500, 0.2500, 0.1750];
params.n_gene=n_gene;
params.P_int=P_k;
PDF =normpdf(x{1}, 270, 25)+normpdf(x{1}, 25, 20); % take the combination of two normal distrubtions
P = PDF/trapz(x{1},PDF); % normalise it, i.e. sum equals to one.
P=P';
%}

u=ones(nt_solve,1)*0;
[du1Lix,P_main]=main_equ(params,u); %solution of the main problem 


 str1 = '#77AC30';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;

str2 = '#A2142F';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;

str3 = '#EDB120';
color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;


str4 = '#0072BD';
color4 = sscanf(str4(2:end),'%2x%2x%2x',[1 3])/255;

%{
figure
hold on
plot(x{1},P_main{end},'-','LineWidth',2.5,'Color',color2)
%plot(x{1},P,'--','LineWidth',1.5)
xlabel('Protein 1')
ylabel('Probability')
hold off

%plot(T,u,'LineWidth',1.5);

%}
%plot screenshot of solution 
%choose time points to plot 
s1=10;
s2=50;
s3=120;
s4=200;
s5=300;

figure 
[X,Y] = meshgrid([s1 s2 s3 s4 s5],x{1});
plot3(X,Y,[ P_main{s1*m} P_main{s2*m} P_main{s3*m} P_main{s4*m} P_main{s5*m}],'Color',color2,LineWidth=2.5)

%plot3(X,Y,[P P P P P],'--','Color',color1,LineWidth=2.5)
xlabel('Time')
ylabel('Protein')
zlabel('Probability')
%legend()

grid on 
view(35,20);

% Save the figure in EPS format
%filename = 'figure_4d.eps';
%print('-depsc', filename);
print('-vector','-dsvg','figure')
%}
