%%%%%%%%% Optimal control for stepwise constant input to solve MPC%%%%%%%%%%%%%%%%%%

function[J,P_main, u_l]=optimal_control_mpc(params,u_)

n_gene=params.n_gene;

Time_mesh=params.time.Time_mesh;
Prot_mesh=params.protein.Prot_mesh;
n_op=params.opt.number_operations;
m=Time_mesh(4);
nt_solve=Time_mesh(3);

T=linspace(Time_mesh(1),Time_mesh(2), Time_mesh(3));

%Target  PDF extraction 
P=params.opt.P_target;


%Inducer bounds
upper_bound=params.opt.Inducer_upbo;
lower_bound=params.opt.Inducer_lobo;



iN=cell(n_gene,1);
x=cell(n_gene,1);
for i=1:n_gene
    iN{i} = Prot_mesh(i,3) + 1;
    x{i} = linspace(Prot_mesh(i,1),Prot_mesh(i,2), iN{i});
end



%Initialization of control input
u_k=u_*0;
u_l=u_;
J=zeros(n_op,1);


%%%%%%%%%%%%%%Gradient calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[du1Lix, P_main]=main_equ(params,u_l); %solution of the main problem 

[P_adj]=adjoint_f(params,P_main,P,u_l); %calculate adjoint_solve solution


%%%%%%%%%%%%Calculate gradient_solve of cost function%%%%%%%%%%%%%%%%%%%%%%%%%
Int_solveegral_1=zeros(nt_solve,1);
Int_solveegral_m=zeros(nt_solve,1);

for i=1:nt_solve
    Int_solveegral_1(i)=trapz(x{1},(du1Lix{i}).*P_adj{i});

end 

for k=1:params.mpc.N_ch
Int_solveegral_m((k-1)*m+1:k*m)=trapz(T((k-1)*m+1:k*m),Int_solveegral_1((k-1)*m+1:k*m));
end


dJ_u1=(Int_solveegral_m);

d1=-dJ_u1;


j=1;

%Initial cost function value 
J(1)=cost_f(P_main,P);




%%%%Running NCG algorithm%%%%%%%
while j<n_op 

%%%%Line search armijo condition%%%%%%
[alpha]=Steplength_Backtracking_m(params,dJ_u1,u_k,u_l,d1,J(j));

%Update control input
u_k=u_l;

u_l=u_l+alpha.*d1;

u_b=u_l>upper_bound;
l_b=u_l<lower_bound;


u_l(u_b)=upper_bound;
u_l(l_b)=lower_bound;


%saving previous solution 
Pk=P_main;

%%%%Evaluation of the new cost function value%%%%%%%%%%%%%%%%%%%%%%%
[du1Lix, P_main]=main_equ(params,u_l); %solution of the main problem 

J(j+1)=cost_f(P_main,P);


%%%%Stopping tests%%%%%%%%
if isnan(J(j+1)) 
     u_l=u_k;
    P_main=Pk;
    disp('Nan Value')
    break 
elseif norm(dJ_u1)==0 
    u_l=u_k;
    P_main=Pk;
    disp('Null gradient')
    break 
elseif J(j+1)<J(1)*10^(-8)
    disp('Tolerance reached')
break
elseif J(j+1)>=J(j)
    u_l=u_k;
    P_main=Pk;
    disp('No decrease')
    break  
end


%%%%%%%Gradient evaluation%%%%%%%%%%%
[P_adj]=adjoint_f(params,P_main,P,u_l);

Int_solveegral_1=zeros(nt_solve,1);

for i=1:nt_solve
    Int_solveegral_1(i)=trapz(x{1},(du1Lix{i}).*P_adj{i});
end 

for k=1:params.mpc.N_ch
Int_solveegral_m((k-1)*m+1:k*m)=trapz(T((k-1)*m+1:k*m),Int_solveegral_1((k-1)*m+1:k*m));
end


g_1k=dJ_u1;


dJ_u1=(Int_solveegral_m);

g_1l=dJ_u1;

yk1=g_1l-g_1k;

L2=trapz(T,d1.*yk1);
L1=(1/L2)*(yk1-2*d1.*trapz(T,yk1.*yk1)/L2);

beta=trapz(T,L1.*g_1l);

d1=-g_1l+beta.*d1;

j=j+1;

       
end 
end 



