
function[J,P_main, u_l]=optimal_control_m(params,u_)

%%%%%%%%% MAIN optimal control PROBLEM%%%%%%%%%%%%%%%%%%
n_gene=params.n_gene;

Time_mesh=params.time.Time_mesh;
Prot_mesh=params.protein.Prot_mesh;
n_op=params.opt.number_operations;

%Extract target 
P=params.opt.P_target;

T=linspace(Time_mesh(1),Time_mesh(2), Time_mesh(3));


iN=cell(n_gene,1);
x=cell(n_gene,1);
for i=1:n_gene
    iN{i} = Prot_mesh(i,3) + 1;
    x{i} = linspace(Prot_mesh(i,1),Prot_mesh(i,2), iN{i});
end



%%%%GRadient_solve descent_solve optimization%%%%%%%%%%%%%%%%%%%%%%%%%%


u_k=u_*0;
u_l=u_;
J=zeros(n_op,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[du1Lix, P_main]=main_equ(params,u_l); %solution of the main problem 

[P_adj]=adjoint_f(params,P_main,P,u_l); %calculate adjoint_solve solution

J(1)=cost_f(P_main,P);

%%%%%%%%%%%%Calculate gradient_solve of cost function%%%%%%%%%%%%%%%%%%%%%%%%%
Int_solveegral_1=zeros(Time_mesh(3),1);

for i=1:Time_mesh(3)
    Int_solveegral_1(i)=trapz(x{1},(du1Lix{i}).*P_adj{i});
end 


dJ_u1=(Int_solveegral_1);

d1=-dJ_u1;

j=1;

upper_bound=params.opt.Inducer_upbo;
lower_bound=params.opt.Inducer_lobo;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve optimization problem 
while j<n_op 

[alpha]=Steplength_Backtracking_m(params,dJ_u1,u_k,u_l,d1,J(j));

%update of inducers
u_k=u_l;

u_l=u_l+alpha.*d1;

u_b=u_l>upper_bound;
l_b=u_l<lower_bound;


u_l(u_b)=upper_bound;
u_l(l_b)=lower_bound;




Pk=P_main;

[du1Lix, P_main]=main_equ(params,u_l); %solution of the main problem 

%%%%%%%%%%%%Calculate gradient_solve of cost function%%%%%%%%%%%%%%%%%%%%%%%%%
J(j+1)=cost_f(P_main,P);

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
elseif J(j+1)<J(1)*params.opt.NCG_tol
    disp('Tolerance reached')
break
elseif J(j+1)>=J(j)
    u_l=u_k;
    P_main=Pk;
    disp('No decrease')
    break 
end


[P_adj]=adjoint_f(params,P_main,P,u_l);

Int_solveegral_1=zeros(Time_mesh(3),1);

for i=1:Time_mesh(3)
    Int_solveegral_1(i)=trapz(x{1},(du1Lix{i}).*P_adj{i});
end 



g_1k=dJ_u1;


dJ_u1=(Int_solveegral_1);

g_1l=dJ_u1;

yk1=g_1l-g_1k;


L2=trapz(T,d1.*yk1);
L1=(1/L2)*(yk1-2*d1.*trapz(T,yk1.*yk1)/L2);

beta=trapz(T,L1.*g_1l);


d1=-g_1l+beta.*d1;



j=j+1;
end 


end 



