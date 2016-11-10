%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function obj=initialize_time_domain_solution(obj,t_step)
obj.Uinf=a2bf(norm(obj.Uinf),obj.state.alpha,obj.state.beta,obj.Ma_corr);
obj.t_step=t_step;
obj.qinf=1/2*obj.state.rho_air*norm(obj.Uinf)^2;
n_wake=obj.settings.wakelength_factor*obj.reference.b_ref/(obj.t_step*norm(obj.Uinf));
obj=obj.initialize_vring_grid(obj.t_step);
obj=obj.initialize_wake_grid(ceil(n_wake),obj.t_step);


x_min=min(obj.grid(1,:));
x_max=max(obj.grid_wake(1,:));
            steps=60;
            obj.x_kin=x_min:obj.Uinf(1)*1/obj.t_step:x_max;
            obj.Uinf_x=zeros(3,steps);
            obj.Uinf_x(1,:)=obj.Uinf(1);
            obj.Uinf_x(2,:)=obj.Uinf(2);
            obj.Uinf_x(3,:)=obj.Uinf(3);


obj.Gamma_wake=zeros(length(obj.panels_wake),1);
obj=obj.compute_influence_coeff_matrix();
obj=obj.determine_boundary_conditions();
obj.Gamma=linsolve(obj.Abb,obj.b');
obj.Gamma_wake=obj.Cbb*obj.Gamma;
obj=obj.shed_wake();

i=1;
while abs((sum(abs(obj.Gamma_prv))-sum(abs(obj.Gamma)))/sum(abs(obj.Gamma)))>10^(-10)
    obj=obj.determine_boundary_conditions();
    obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
    bx=obj.b-(obj.Abw*obj.Gamma_wake)';
    obj.Gamma_prv=obj.Gamma;
    obj.Gamma=linsolve(obj.Abb,bx');
    obj=obj.shed_wake();
    i=i+1;
    % abort for i>limit ?
end
obj=obj.f_postprocess;
mkdir(['results/' obj.case_name],'free_flying_mean_axis');
end
