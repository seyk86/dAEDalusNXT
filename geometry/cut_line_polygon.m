function [new_polygon ] =cut_line_polygon(master_quad,line)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
k=0;

new_polygon=zeros(3,1);
newpoints=zeros(3,1);

point_face=0;
point_dist=0;


diff_1=line(:,1)-line(:,2);


[nono,comp_idx_2]=max(diff_1);
diff_1(comp_idx_2)=0;
[nono,comp_idx_1]=max(diff_1);
if comp_idx_1==comp_idx_2
    if comp_idx_1==3
        comp_idx_2=1;
    else
        comp_idx_2=comp_idx_1+1;
    end
end

n=size(master_quad,2);
for i=1:n
    O1=master_quad(:,i);
    if i<n
        O2=master_quad(:,i+1);
    else
        O2=master_quad(:,1);
    end
    
    ro=O2-O1;
    normro=norm(ro);
    
    I1=line(:,1);
    I2=line(:,2);
    ri=I2-I1;
    normri=norm(ri);
    if ~((norm(cross(ro,ri))< eps) || (abs(ro(1)/normro-ri(1)/normri)<10*eps &&  abs(ro(2)/normro-ri(2)/normri)<10*eps && abs(ro(3)/normro-ri(3)/normri)<10*eps))
                    lambda=((O1(comp_idx_2)-I1(comp_idx_2)+I1(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1)-O1(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1))/(-ri(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1)+ri(comp_idx_2)));
                    xi=(I1(comp_idx_1)-O1(comp_idx_1)+lambda*ri(comp_idx_1))/ro(comp_idx_1);
        %U1=I1+lambda*ri
        U1=O1+xi*ro;
        % check if point interpolated ok

        % check if point is in between P1 and P2
        %lambda=(U1(1)-I1(1))/(I2(1)-I1(1));
        %errU1_online=norm(U1-I1-lambda*(I2-I1));
        if lambda<1.0-100*eps && lambda>0.0+100*eps && xi<1.0-100*eps && xi>0.0+100*eps
                    errU1=norm(U1-I1-lambda*ri);
%         if errU1>1E-10
%             disp('geometry error');
%         end
            
            k=k+1;
            newpoints(:,k)=U1;
            point_face(k)=1;
            point_dist(k)=lambda;
        end
    end
end
if k>0
    new_polygon=zeros(3,k);
    ll=1;
    
    tosort=[point_face;point_dist;1:k];
    sorted=sortrows(tosort')';
    for j=1:2
        if ll>=2
           if ~((abs(new_polygon(1,ll-1)-line(1,j))<0.001) && (abs(new_polygon(2,ll-1)-line(2,j))<0.001) &&(abs(new_polygon(3,ll-1)-line(3,j))<0.001))
               new_polygon(:,ll)=line(:,j);
               ll=ll+1;
           end 
        else
            new_polygon(:,ll)=line(:,j);
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
    new_polygon=line;
end




fact=0.0;
end



