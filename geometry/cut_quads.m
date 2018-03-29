function [new_polygon ] =cut_quads(master_quad,panel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
k=0;

new_polygon=zeros(3,1);
newpoints=zeros(3,1);

point_face=0;
point_dist=0;

for i=1:4
    O1=master_quad(:,i);
    if i<4
        O2=master_quad(:,i+1);
    else
        O2=master_quad(:,1);
    end
    
    ro=O2-O1;
    normro=norm(ro);
    
    
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
                    xi=(I1(1)-O1(1)+lambda*ri(1))/ro(1);
                    U1=I1+lambda*ri;
                    % check if point interpolated ok
                    errU1=norm(U1-I1-lambda*ri);
                    if errU1>1E-10
                        disp('geometry error');
                    end
                    % check if point is in between P1 and P2
                    lambda=(U1(1)-I1(1))/(I2(1)-I1(1));
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
                    xi=(I1(1)-O1(1)+lambda*ri(1))/ro(1);
                    U1=I1+lambda*ri;
                    % check if point interpolated ok
                    errU1=norm(U1-I1-lambda*ri);
                    if errU1>1E-10
                        disp('geometry error');
                    end
                    % check if point is in between P1 and P2
                    lambda=(U1(1)-I1(1))/(I2(1)-I1(1));
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
end
if k>0
    new_polygon=zeros(3,k);
    ll=1;
    
    tosort=[point_face;point_dist;1:k];
    sorted=sortrows(tosort')';
    for j=1:4
        if ll>=2
           if ~((abs(new_polygon(1,ll-1)-panel(1,j))<0.001) && (abs(new_polygon(2,ll-1)-panel(2,j))<0.001) &&(abs(new_polygon(3,ll-1)-panel(3,j))<0.001))
               new_polygon(:,ll)=panel(:,j);
               ll=ll+1;
           end 
        else
            new_polygon(:,ll)=panel(:,j);
            ll=ll+1;
        end
        for k=1:length(point_face)
            if sorted(1,k)==j
               if ~((abs(new_polygon(1,ll-1)-newpoints(1,sorted(3,k)))<0.001) && (abs(new_polygon(2,ll-1)-newpoints(2,sorted(3,k)))<0.001) &&(abs(new_polygon(3,ll-1)-newpoints(3,sorted(3,k)))<0.001))
                    new_polygon(:,ll)=newpoints(:,sorted(3,k));
                    ll=ll+1;
                end
            end
        end 
    end   
    
else
    new_polygon=panel;
end




fact=0.0;
end



