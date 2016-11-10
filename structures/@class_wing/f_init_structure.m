%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%

% Author:           Klaus Seywald
%                   klaus.seywald@mytum.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[wing_struct]=f_init_structure(wing_struct,wing_geo,varargin)
is_wetted=[];
if nargin==3;
    is_wetted=varargin{1};
end

node_coords_x=[];
node_coords_y=[];
node_coords_z=[];
is_fueled_per_element=[];
i=1;
if length(wing_geo.wing_segments)>1
    for i=1:length(wing_geo.wing_segments)-1
        if ~isempty(is_wetted)
            is_fueled_per_element=[is_fueled_per_element ones(1,length(squeeze(wing_geo.wing_segments(i).wingbox_coords(1,1:end-1,1))))*is_wetted(i)];
        end
        node_coords_x=[node_coords_x 0.5*(wing_geo.wing_segments(i).wingbox_coords(1,1:end-1,1)+wing_geo.wing_segments(i).wingbox_coords(1,1:end-1,2))];
        node_coords_y=[node_coords_y 0.5*(wing_geo.wing_segments(i).wingbox_coords(2,1:end-1,1)+wing_geo.wing_segments(i).wingbox_coords(2,1:end-1,2))];
        node_coords_z=[node_coords_z 0.5*(wing_geo.wing_segments(i).wingbox_coords(3,1:end-1,1)+wing_geo.wing_segments(i).wingbox_coords(3,1:end-1,2))];
    end
    i=i+1;
    if ~isempty(is_wetted)
        is_fueled_per_element=[is_fueled_per_element ones(1,length(squeeze(wing_geo.wing_segments(i).wingbox_coords(1,1:end-1,1))))*is_wetted(i)];
    end
else
    if ~isempty(is_wetted)
        is_fueled_per_element=[is_fueled_per_element ones(1,length(squeeze(wing_geo.wing_segments(i).wingbox_coords(1,1:end-1,1))))*is_wetted(i)];
    end

end
node_coords_x=[node_coords_x 0.5*(wing_geo.wing_segments(end).wingbox_coords(1,:,1)+wing_geo.wing_segments(end).wingbox_coords(1,:,2))];
node_coords_y=[node_coords_y 0.5*(wing_geo.wing_segments(end).wingbox_coords(2,:,1)+wing_geo.wing_segments(end).wingbox_coords(2,:,2))];
node_coords_z=[node_coords_z 0.5*(wing_geo.wing_segments(end).wingbox_coords(3,:,1)+wing_geo.wing_segments(end).wingbox_coords(3,:,2))];

if wing_geo.symmetric==1
    if ~isempty(is_wetted)
        is_fueled_per_element=[is_fueled_per_element(end:-1:1) is_fueled_per_element(1:end)];
    end
    node_coords_x=[node_coords_x(end:-1:1) node_coords_x(2:end)];
    node_coords_y=[-node_coords_y(end:-1:1) node_coords_y(2:end)];
    node_coords_z=[node_coords_z(end:-1:1)  node_coords_z(2:end)];
end


wing_struct.node_coords(:,1)=node_coords_x;
wing_struct.node_coords(:,2)=node_coords_y;
wing_struct.node_coords(:,3)=node_coords_z;

wing_struct.wing_frontview_length=0;
for i=1:length(wing_geo.wing_segments)
    wing_struct.wing_frontview_length=wing_struct.wing_frontview_length+wing_geo.wing_segments(i).b;
end

wing_struct.r_Ref(1)=wing_geo.wingbox_coords(1,1,1);
wing_struct.r_Ref(2)=wing_geo.wingbox_coords(2,1,1);
wing_struct.r_Ref(3)=wing_geo.wingbox_coords(3,1,1);

wing_struct.is_sym=wing_geo.symmetric;
wing_struct.isExternalFEM = wing_geo.isExternalFEM;

%% initialize wing segment structure
%crosssection=class_crosssection_wingbox;

k=1;
V_wingbox=0;
S_wing=0;
for i=1:length(wing_geo.wing_segments)
    span_loc_grid=0:1/(size(wing_geo.wing_segments(i).wingbox_coords,2)-1):1;
    for j=1:size(wing_geo.wing_segments(i).wingbox_coords,2)-1
        if wing_geo.isExternalFEM==0
            c=wing_geo.wing_segments(i).c_r*(1-span_loc_grid(j+1)*0.5-span_loc_grid(j)*0.5)+wing_geo.wing_segments(i).c_t*(span_loc_grid(j+1)*0.5+span_loc_grid(j)*0.5);
            h=0.25*(wing_geo.wing_segments(i).wingbox_height(j,1)+wing_geo.wing_segments(i).wingbox_height(j,2)+wing_geo.wing_segments(i).wingbox_height(j+1,1)+wing_geo.wing_segments(i).wingbox_height(j+1,2));
        end
      %  h=0.013;
        w=0.5*(norm(wing_geo.wing_segments(i).wingbox_coords(:,j,1)-wing_geo.wing_segments(i).wingbox_coords(:,j,2))+norm(wing_geo.wing_segments(i).wingbox_coords(:,j+1,1)-wing_geo.wing_segments(i).wingbox_coords(:,j+1,2)));
 
        
        elem_vec=0.5*(wing_geo.wing_segments(i).wingbox_coords(:,j+1,1)+wing_geo.wing_segments(i).wingbox_coords(:,j+1,2))...
            -0.5*(wing_geo.wing_segments(i).wingbox_coords(:,j,1)+wing_geo.wing_segments(i).wingbox_coords(:,j,2));
        
        
        le=sqrt(elem_vec(1)*elem_vec(1)+elem_vec(2)*elem_vec(2)+elem_vec(3)*elem_vec(3));
        phi=-asin(elem_vec(1)/le);
        if abs(elem_vec(2))<1E-9
            if elem_vec(3)>0
                nu=90*pi/180;
            elseif elem_vec(3)<0
                nu=270*pi/180;
            end
        else
            if elem_vec(2)>0
                nu=atan(elem_vec(3)/elem_vec(2));
            elseif elem_vec(2)<0
                nu=atan(elem_vec(3)/elem_vec(2))+pi;
            end
        end
        if ((elem_vec(2)==0) && (elem_vec(3)==0))
           nu=0;
           twist=0;
%            if wing_geo.wing_segments(i).wingbox_coords(2,j,1)<0
%                phi=-phi;
%            end
            epsilon(k)=0;
        else
            twist=(wing_geo.wing_segments(i).Theta_r*(1-((span_loc_grid(j+1)*0.5+span_loc_grid(j)*0.5)))+wing_geo.wing_segments(i).Theta_t*(span_loc_grid(j+1)*0.5+span_loc_grid(j)*0.5))*pi/180;
            epsilon(k)=(wing_geo.wing_segments(i).Theta_r*(1-span_loc_grid(j))+wing_geo.wing_segments(i).Theta_t*(span_loc_grid(j)))*pi/180;
        end
        %%TODO fix:
        if abs(phi+pi/2)<1E-6
           nu=pi/2;
           %epsilon(k)=pi/2;
        end
        if abs(phi-pi/2)<1E-6
           nu=pi/2;
          % epsilon(k)=pi/2;
        end
        
        
        nnu(k)=nu;
        nphi(k)=phi;
        %%
        if wing_geo.isExternalFEM==0
            V_wingbox=V_wingbox+h*w*le;
            S_wing=S_wing+0.8*(2*h+2*c)*le;
                
            wing_struct.beamelement(k).crosssection=wing_struct.beamelement(k).crosssection.setGeometry(c,h,w);%*cos(phi));
        end
        if wing_geo.symmetric==1
            wing_struct.beamelement(k)=wing_struct.beamelement(k).setElementGeometry(le,-phi,-nu,twist);
        else
            wing_struct.beamelement(k)=wing_struct.beamelement(k).setElementGeometry(le,phi,nu,twist);
        end
        
        k=k+1;
        %% todo calculate more accurate
       
    end
    if wing_geo.isExternalFEM==0
        wing_struct.Awing=wing_struct.Awing+wing_geo.wing_segments(i).S;
        wing_struct.V_wingbox=V_wingbox;
    end
end


if wing_geo.symmetric==1
    for i=1:length(wing_geo.wing_segments)
        span_loc_grid=0:1/(size(wing_geo.wing_segments(i).wingbox_coords,2)-1):1;
        for j=1:size(wing_geo.wing_segments(i).wingbox_coords,2)-1
            if wing_geo.isExternalFEM==0
                c=wing_geo.wing_segments(i).c_r*(1-span_loc_grid(j+1)*0.5-span_loc_grid(j)*0.5)+wing_geo.wing_segments(i).c_t*(span_loc_grid(j+1)*0.5+span_loc_grid(j)*0.5);
                h=0.25*(wing_geo.wing_segments(i).wingbox_height(j,1)+wing_geo.wing_segments(i).wingbox_height(j,2)+wing_geo.wing_segments(i).wingbox_height(j+1,1)+wing_geo.wing_segments(i).wingbox_height(j+1,2));
        %    h=0.013;
            end
            w=0.5*(norm(wing_geo.wing_segments(i).wingbox_coords(:,j,1)-wing_geo.wing_segments(i).wingbox_coords(:,j,2))+norm(wing_geo.wing_segments(i).wingbox_coords(:,j+1,1)-wing_geo.wing_segments(i).wingbox_coords(:,j+1,2)));
            wing_struct.beamelement(k)=wing_struct.beamelement(k-1);
            
            elem_vec=0.5*(wing_geo.wing_segments(i).wingbox_coords(:,j+1,1)+wing_geo.wing_segments(i).wingbox_coords(:,j+1,2))...
                -0.5*(wing_geo.wing_segments(i).wingbox_coords(:,j,1)+wing_geo.wing_segments(i).wingbox_coords(:,j,2));
 
            le=sqrt(elem_vec(1)*elem_vec(1)+elem_vec(2)*elem_vec(2)+elem_vec(3)*elem_vec(3));
            phi=-asin(elem_vec(1)/le);
            if abs(elem_vec(2))<1E-9
                if elem_vec(3)>0
                    nu=90*pi/180;
                elseif elem_vec(3)<0
                    nu=270*pi/180;
                end
            else
                if elem_vec(2)>0
                    nu=atan(elem_vec(3)/elem_vec(2));
                elseif elem_vec(2)<0
                    nu=atan(elem_vec(3)/elem_vec(2))+pi;
                end
            end
            
            twist=(wing_geo.wing_segments(i).Theta_r*(1-((span_loc_grid(j+1)*0.5+span_loc_grid(j)*0.5)))+wing_geo.wing_segments(i).Theta_t*(span_loc_grid(j+1)*0.5+span_loc_grid(j)*0.5))*pi/180;
            %%TODO fix:      
            %%
            if wing_geo.isExternalFEM==0
                V_wingbox=V_wingbox+h*w*le;
                        
                wing_struct.beamelement(k).crosssection=wing_struct.beamelement(k).crosssection.setGeometry(c,h,w);            
            end
            wing_struct.beamelement(k)=wing_struct.beamelement(k).setElementGeometry(le,phi,nu,twist);
            k=k+1;
            %% todo calculate more accurate
          
        end
          wing_struct.Awing=wing_struct.Awing+wing_geo.wing_segments(i).S;
          wing_struct.V_wingbox=V_wingbox;
    end
%     for i=1:wing_struct.nel/2
%        zwsp(i)=wing_struct.beamelement(i);
%     end
%     for i=1:wing_struct.nel/2
%        wing_struct.beamelement(i)=zwsp(wing_struct.nel/2-i+1);
%     end
%aircraft_structure.beam(2).epsilon=[aircraft_structure.beam(2).epsilon(end/2-1:-1:1) 0 0 aircraft_structure.beam(2).epsilon(1:end/2)];
%aircraft_structure.beam(2).nu=[-aircraft_structure.beam(2).nu(end/2-1:-1:1) 0 0 aircraft_structure.beam(2).nu(1:end/2)];
%aircraft_structure.beam(2).phi=[-aircraft_structure.beam(2).phi(end/2-1:-1:1) 0 0 aircraft_structure.beam(2).phi(1:end/2)];

    wing_struct.beamelement(1:wing_struct.nel/2)=wing_struct.beamelement(wing_struct.nel/2:-1:1);
end

if ((elem_vec(2)==0) && (elem_vec(3)==0))
    nu=0;
    twist=0;
    epsilon(k)=0;
else
    epsilon(k)=wing_geo.wing_segments(end).Theta_t*pi/180;

end
    nnu(k)=nu;
nphi(k)=phi;
dist_c4_sc=zeros(length(wing_struct.beamelement)+1,1);
nn=1;
for i=1:1:length(wing_struct.beamelement)+1
    if wing_geo.symmetric==1
        if i<=length(wing_geo.wingbox_c4(1,:))
            dist_c4_sc(i)=norm(wing_geo.wingbox_c4(:,end+1-i)'-wing_struct.node_coords(end+1-i,:));
            nn=i;
        else
            wing_geo.wingbox_c4(2,:)=-wing_geo.wingbox_c4(2,:);
            dist_c4_sc(i)=norm(wing_geo.wingbox_c4(:,i-nn+1)'-wing_struct.node_coords(end+1-i,:));
            wing_geo.wingbox_c4(2,:)=-wing_geo.wingbox_c4(2,:);
            %dist_c4_sc(i)=norm(wing_geo.wingbox_c4(:,i-nn)'-wing_struct.node_coords(i-nn,:));
        end
    else
        dist_c4_sc(i)=norm(wing_geo.wingbox_c4(:,i)'-wing_struct.node_coords(i,:));
    end 
end

wing_struct.dist_c4_sc=dist_c4_sc;
wing_struct.epsilon=epsilon;
wing_struct.nu=nnu;
wing_struct.phi=nphi;


if~isempty(is_fueled_per_element)
    for i=1:length(wing_struct.beamelement)
        wing_struct.beamelement(i).is_fueled=is_fueled_per_element(i);
    end
end

end % function
