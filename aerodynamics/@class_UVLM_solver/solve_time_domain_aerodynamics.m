%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function obj=solve_time_domain_aerodynamics(obj,aircraft,x_body,V,alpha,beta,pqr,rho_air,itr,filename)
obj.Uinf=a2bf(norm(V),alpha,beta,obj.Ma_corr);
obj.qinf=1/2*rho_air*norm(obj.Uinf)^2;

obj.Uinf_x(1,:)=obj.Uinf(1);
obj.Uinf_x(2,:)=obj.Uinf(2);
obj.Uinf_x(3,:)=obj.Uinf(3);

prv_colloc=obj.colloc;
if isempty(aircraft.grid_deflected)
    obj.grid(1,:)=aircraft.grid(1,:)/obj.Ma_corr;
    obj.grid(2,:)=aircraft.grid(2,:);
    obj.grid(3,:)=aircraft.grid(3,:);
else
    obj.grid(1,:)=aircraft.grid_deflected(1,:)/obj.Ma_corr;
    obj.grid(2,:)=aircraft.grid_deflected(2,:);
    obj.grid(3,:)=aircraft.grid_deflected(3,:);
end
obj=obj.update_grid();
obj=obj.determine_boundary_conditions_gust_deflection_rotation(pqr(1),pqr(2),pqr(3),prv_colloc,obj.t_step);
obj.Gamma_wake=obj.Cbb*obj.Gamma+[zeros(obj.panel_len,1); obj.Gamma_wake(1:end-obj.panel_len)];
bx=obj.b-(obj.Abw*obj.Gamma_wake)';
obj.Gamma_prv=obj.Gamma;
obj.Gamma=linsolve(obj.Abb,bx');
obj=obj.f_postprocess;

% if gust==1
%     obj=obj.shed_wake_gust();
% else
    obj=obj.shed_wake();
% end


%obj=obj.perform_wake_rollup;

%obj=obj.compute_influence_coeff_matrix();

if obj.settings.movie==1
    obj.write_tecplot_free_flying_wake([filename num2str(itr)],[x_body(1:3)'; x_body(4:6)'],1);
end

end
