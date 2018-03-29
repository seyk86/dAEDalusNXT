classdef class_nacellesegment<class_circularsegment
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        S_wet;
    end
    
    methods
        function obj = class_nacellesegment(pos,a,b,l,np,taper_f,taper_r,varargin)
            swn=[];
            sw=[];
            if nargin==9
                swn=varargin{1};
                sw=varargin{2};
            end
            obj=obj@class_circularsegment(pos,a,b,l,np,taper_f,taper_r,swn,sw); 
        end
        
        function obj= compute_wetted_area(obj)
            A=0;
            
            for i=1:2:length(obj.xyz)-3
  
             v1=obj.xyz(:,i+1)-obj.xyz(:,i);
             v2=obj.xyz(:,i+2)-obj.xyz(:,i);
             A1=norm(cross(v1,v2))/2;

             v3=obj.xyz(:,i+1)-obj.xyz(:,i+3);
             v4=obj.xyz(:,i+2)-obj.xyz(:,i+3);
             A2=norm(cross(v3,v4))/2;
             A=A+A1+A2;
            end
            obj.S_wet=A;
        end
        
    end
    
end

