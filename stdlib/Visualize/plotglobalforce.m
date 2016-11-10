function [ output_args ] = plotglobalforce(handle,FM,scale,node_coords,direction,color,linecolor,flag_plotbeam)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
            
            hold on
            grid on
            if(flag_plotbeam)  
                plot3(handle,node_coords(:,1),node_coords(:,2),node_coords(:,3),'-k','LineWidth',2)
            end
            %if(direction=='z')
             %   plot3(handle,node_coords(:,1),node_coords(:,2),node_coords(:,3)+FM./scale,'r-','LineWidth',1)
            %end
            
            for i=1:1:(length(FM)-1)
                
                x=node_coords(i,1);
                y=node_coords(i,2);
                z=node_coords(i,3);
                
                xp=node_coords(i+1,1);
                yp=node_coords(i+1,2);
                zp=node_coords(i+1,3);
                
                Mx=FM(i);
                Mxp=FM(i+1);
                
                if(direction=='x')
                    pvec=[0,0,-1];
                end
                
                if(direction=='y')
                    pvec=[0,1,0];
                end
                
                if(direction=='z')
                   pvec=[1,0,0];
                end
                
                vec1=[0,yp-y,zp-z];  
                vec2=cross(pvec,vec1);

                %vec2=cross(vec1,normv);
                vec2=vec2/norm(vec2);
                   
                vec3=vec2*Mx/scale;
                vec4=vec2*Mxp/scale; 
                    
                X=[x,xp,xp+vec4(1),x+vec3(1),x+vec3(1),x+vec3(1),x+vec3(1),x+vec3(1)];
                Y=[y,yp,yp+vec4(2),y+vec3(2),y+vec3(2),y+vec3(2),y+vec3(2),y+vec3(2)];
                Z=[z,zp,zp+vec4(3),z+vec3(3),z+vec3(3),z+vec3(3),z+vec3(3),z+vec3(3)];     

                line([xp+vec4(1) x+vec3(1)],[yp+vec4(2),y+vec3(2)],[zp+vec4(3),z+vec3(3)],'Color','k'); 
                line([x x+vec3(1)],[y y+vec3(2)],[z z+vec3(3)],'Color',linecolor); 
                line([xp xp+vec4(1)],[yp yp+vec4(2)],[zp zp+vec4(3)],'Color',linecolor); 
                %line([xp+vec4(1) x+vec3(1)],[yp+vec4(2),y+vec3(2)],[zp+vec4(3),z+vec3(3)],'Color','k'); 
                hdl=fill3(X,Y,Z,[color color color color color 0 1 color]);
                %set transparency
                alpha(hdl,0.2)
               % set(hdl,'EdgeColor','r'); 
                set(hdl,'LineStyle','none');     
            end
            axis equal
end

