classdef class_nacelle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties
        nacelle_segments;
        
        len;
        diam_in;
        diam_max;
        diam_out;
        CD_f;
        D_f;
        S_wet;
        
        grid;
        panels;
        te_idx;
    end
    
    methods
        function obj=class_nacelle(pos,a,b,l,np)
            if length(a)==length(b) &&( length(b)==(length(l)+1))
                obj.len=0;
                obj.diam_in=a(1)+b(1);
                obj.diam_out=a(end)+b(end);
                obj.diam_max=max(a+b);
                for i=1:length(a)-1
                    newsegment=class_nacellesegment(pos,a(i),b(i),l(i),np,1,a(i+1)/a(i));
                    newsegment.n_seg=3;
                    obj.nacelle_segments=[obj.nacelle_segments newsegment];
                    pos=pos+[l(i) 0 0];
                    obj.len=obj.len+l(i);
                end
            end
            obj=obj.compute_wetted_area();
        end
        
        function obj=compute_wetted_area(obj)
            obj.S_wet=0;
            for i=1:length(obj.nacelle_segments)
                obj.nacelle_segments(i)=obj.nacelle_segments(i).compute_wetted_area();
                obj.S_wet=obj.S_wet+obj.nacelle_segments(i).S_wet;
            end
        end
        
        function obj=plot_nacelle(obj)
            for i=1:length(obj.nacelle_segments)
                obj.nacelle_segments(i).plot_circularsegment();
            end
            axis equal
        end
        
        function obj=compute_friction_drag(obj,state,S_ref)
            Re_l=state.rho_air*norm(state.V_inf)*obj.len/state.mu;
            % laminar or turbulent
            if Re_l>5E5
                CF=0.455/log10(Re_l)^2.58;
            else
                CF=1.328/sqrt(Re_l);
            end
            % form factor for nacelles: http://adg.stanford.edu/aa241/drag/formfactor.html
            Amax=obj.diam_max*pi;
            Aexit=obj.diam_out*pi;
            Ainflow=obj.diam_in*pi;
            LD=(obj.len+obj.diam_in)/sqrt(4/pi*(Amax-(Aexit+Ainflow)/2));
            D=(1-(1-state.Ma^2)/LD^2)^0.5;
            a=2*(1-state.Ma^2)/(LD^2*D^3)*(atanh(D)-D);
            du_maxU0=a/(2-a)/(1-state.Ma^2)^0.5;
            C=2.3;
            form_factor=(1+C*du_maxU0)^2;
            
            % nacelle base drag: http://adg.stanford.edu/aa241/drag/nacbasedrag.html
            
            CD_nacellebase=0.5/12*pi*obj.diam_out*0.2/S_ref;
            obj.CD_f=form_factor*CF*obj.S_wet/S_ref+CD_nacellebase;
            obj.D_f=form_factor*1/2*state.rho_air*norm(state.V_inf)^2*CF*obj.S_wet;
        end
        
        function obj=compute_grid(obj)
            grid=[];
            te_idx=[];
            panels=[];
            len=[];
            
            for i=1:length(obj.nacelle_segments)
                obj.nacelle_segments(i)=obj.nacelle_segments(i).compute_grid();
                len=length(grid);
                grid=[grid obj.nacelle_segments(i).grid];
                te_idx=[te_idx obj.nacelle_segments(i).te_idx+len];
                panels=[panels,obj.nacelle_segments(i).panels+len];
            end
            obj.grid=grid;
            obj.te_idx=te_idx;
            obj.panels=panels;
            
        end
        
        function obj=plot_grid(obj)
            hold on
            for i=1:length(obj.panels)
                handle= fill3(obj.grid(1,obj.panels(:,i)), obj.grid(2,obj.panels(:,i)),obj.grid(3,obj.panels(:,i)),'b');
                alpha(handle,0.4);
            end
        end
    end
end
    
