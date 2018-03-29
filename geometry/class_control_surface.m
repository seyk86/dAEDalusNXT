classdef class_control_surface
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> control surface name
        name;
        %>  1 leading then edge devic, 0 trailing edge device
        pos = 0; 
        %> tapered or constant
        is_tapered = 0;
        %> number of hinge lines
        n_hinge;
        %> chord
        hinge;
        %> deflection
        delta;
        
        delta_l_r;
        %> is symmetric
        is_sym;
        %> deflection
        is_sym_defl;
    end
    
    methods
        function obj=class_control_surface(name,pos,hinge,varargin)
            obj.name=name;
            obj.pos=pos;
            obj.hinge=str2num(hinge);
            obj.n_hinge=length(obj.hinge);
            obj.delta=zeros(1,obj.n_hinge);
            if nargin==4
                obj.is_sym=1;
                obj.is_sym_defl=varargin{1};
                
                if obj.is_sym_defl==0
                    obj.delta_l_r=[0 0];
                end
            end
            if nargin==5
                obj.is_sym=1;
                obj.is_sym_defl=varargin{1};
                obj.is_tapered=varargin{2};
                if obj.is_sym_defl==0
                    obj.delta_l_r=[0 0];
                end
            end
            
        end
    end
    
end

