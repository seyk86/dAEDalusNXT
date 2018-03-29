function [new_polygon ] =cut_line_quad(master_line,panel)
% load debug_vars
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% figure;
% hold on
% scatter3(panel(1,:),panel(2,:),panel(3,:), 's', 'filled')
% h=fill3(panel(1,:),panel(2,:),panel(3,:),'g')
% text(panel(1,1),panel(2,1),panel(3,1),'P1');
% text(panel(1,2),panel(2,2),panel(3,2),'P2');
% text(panel(1,3),panel(2,3),panel(3,3),'P3');
% text(panel(1,4),panel(2,4),panel(3,4),'P4');
% set(h,'facealpha',0.1)
% line(master_line(1,:),master_line(2,:),master_line(3,:), 'Linewidth', 3,'Color', 'r')

k=0;

new_polygon=zeros(3,1);
newpoints=zeros(3,1);

point_face=0;
point_dist=0;
% normal of panel
diagPanel1=panel(:,3)-panel(:,1);
diagPanel2=panel(:,4)-panel(:,2);
normalPanel=cross(diagPanel1,diagPanel2);
normalPanel=normalPanel./norm(normalPanel);


O1=master_line(:,1);
O2=master_line(:,2);
%vector of cutting line
ro=O2-O1;
normro=norm(ro);
lineVector=ro./normro;
%normal of "cutting" plane containing master_line (lies within it) and normal to panel
normalCuttingPlane=cross(normalPanel,lineVector);


for j=1:4
    if j<4
        %   if (in(j)==1 && in(j+1)==1)%||(in(j)==0 && in(j+1)==0)
        
        %else
        % check if on line
        I1=panel(:,j);
        I2=panel(:,j+1);
        ri=I2-I1;
        
        normri=norm(ri);
        if ~((norm(cross(ro,ri))< eps) || (abs(ro(1)/normro-ri(1)/normri)<10*eps &&  abs(ro(2)/normro-ri(2)/normri)<10*eps && abs(ro(3)/normro-ri(3)/normri)<10*eps))
            lambda=((O1(2)-I1(2)+I1(1)*ro(2)/ro(1)-O1(1)*ro(2)/ro(1))/(-ri(1)*ro(2)/ro(1)+ri(2)));
            %simon: better way for finding intersection of line and line is
            %by using a plane containing one of the lines normal to the
            %other line
            lambda=dot((O1-I1),normalCuttingPlane)/dot(ri,normalCuttingPlane);
            
            xi=(I1(1)-O1(1)+lambda*ri(1))/ro(1); %defines the lambda for the line O1 to O2 -> not really necessary?
            U1=I1+lambda*ri;
            % check if point interpolated ok
            errU1=norm(U1-I1-lambda*ri);
            if errU1>1E-10
                disp('geometry error');
            end
            % check if point is in between P1 and P2
            %lambda=(U1(1)-I1(1))/(I2(1)-I1(1)); %SIMON: THIS DOES NOT
            %WORK!!!!!
            %better: use the lambda as calculated above -> it is already
            %correct when between 0 and 1
            
            errU1_online=norm(U1-I1-lambda*(I2-I1));
            if lambda<1.0+100*eps && lambda>0.0-100*eps && xi<1.0+100*eps && xi>0.0-100*eps
                k=k+1;
                newpoints(:,k)=U1;
                point_face(k)=j;
                point_dist(k)=lambda;
            end
        end
        %  end
    else
        %  if (in(j)==1 && in(1)==1) %||(in(j)==0 && in(1)==0)
        
        %  else
        I1=panel(:,j);
        I2=panel(:,1);
        ri=I2-I1;
        normri=norm(ri);
        if ~((norm(cross(ro,ri))< eps) || (abs(ro(1)/normro-ri(1)/normri)<10*eps &&  abs(ro(2)/normro-ri(2)/normri)<10*eps && abs(ro(3)/normro-ri(3)/normri)<10*eps))
            lambda=((O1(2)-I1(2)+I1(1)*ro(2)/ro(1)-O1(1)*ro(2)/ro(1))/(-ri(1)*ro(2)/ro(1)+ri(2)));
            lambda=dot((O1-I1),normalCuttingPlane)/dot(ri,normalCuttingPlane); %see above for comments
            xi=(I1(1)-O1(1)+lambda*ri(1))/ro(1);
            U1=I1+lambda*ri;
            % check if point interpolated ok
            errU1=norm(U1-I1-lambda*ri);
            if errU1>1E-10
                disp('geometry error');
            end
            % check if point is in between P1 and P2
            %lambda=(U1(1)-I1(1))/(I2(1)-I1(1));
            errU1_online=norm(U1-I1-lambda*(I2-I1));
            if lambda<1.0+100*eps && lambda>0.0-100*eps && xi<1.0+100*eps && xi>0.0-100*eps
                k=k+1;
                newpoints(:,k)=U1;
                point_face(k)=j;
                point_dist(k)=lambda;
            end
        end
        % end
    end
end
if k>0
    new_polygon=zeros(3,k);
    ll=1;
    
    tosort=[point_face;point_dist;1:k];
    sorted=sortrows(tosort')';
    for j=1:4
        for k=1:length(point_face)
            if sorted(1,k)==j
                if ll>=2
                    if ~((abs(new_polygon(1,ll-1)-newpoints(1,sorted(3,k)))<100*eps) && (abs(new_polygon(2,ll-1)-newpoints(2,sorted(3,k)))<100*eps) &&(abs(new_polygon(3,ll-1)-newpoints(3,sorted(3,k)))<100*eps))
                        new_polygon(:,ll)=newpoints(:,sorted(3,k));
                        ll=ll+1;
                    end
                else
                    new_polygon(:,ll)=newpoints(:,sorted(3,k));
                    ll=ll+1;
                end
            end
        end
    end
    
else
    new_polygon=panel;
end
%scatter3(new_polygon(1,:),new_polygon(2,:),new_polygon(3,:),'sb', 'filled');


end



