function [newpoints,point_face ] =cut_polygons(master_quad,panel)
%UNTITLED2 Summary of this function goes here
%load('debug_cut_polygon2.mat')
%   Detailed explanation goes here
k=0;
debug=0;
new_polygon=zeros(3,1);
newpoints=zeros(3,1);

point_face=0;
point_dist=0;

diff_1=panel(:,3)-panel(:,1);   %diag 1
diff_2=panel(:,2)-panel(:,4);   %diag 2
diff_1=0.5*diff_1+0.5*diff_2;   %represents leading edge of panel (P2-P1)


[nono,comp_idx_2]=max(diff_1); %first main direction
diff_1(comp_idx_2)=0;
[nono,comp_idx_1]=max(diff_1); %second main direction
if comp_idx_1==comp_idx_2
    if comp_idx_1==3
        comp_idx_2=1;
    else
        comp_idx_2=comp_idx_1+1;
    end
end

n=size(master_quad,2);
for i=1:n   %for every face of the block
    O1=master_quad(:,i);   
    if i<n
        O2=master_quad(:,i+1);
        if i<n-1
            O3=master_quad(:,i+2);
        else
            O3=master_quad(:,1);
        end
    else
        O2=master_quad(:,1);
        O3=master_quad(:,2);
    end
    
    ro=O2-O1;    %face direction of the block
    normro=norm(ro);
    seg_norm=cross((O2-O1),(O3-O1)); %normal vector of the block
    cut_norm=cross(seg_norm,(O2-O1)); % normal vector of cutting plane normal to block and containing current face of block
    
    for j=1:4  %for every face of the panel
        if j<4
         %   if (in(j)==1 && in(j+1)==1)%||(in(j)==0 && in(j+1)==0)
                
            %else
                % check if on line
                I1=panel(:,j);
                I2=panel(:,j+1);
                ri=I2-I1;
                normri=norm(ri);
                if ~((norm(cross(ro,ri))< eps) || (abs(ro(1)/normro-ri(1)/normri)<10*eps &&  abs(ro(2)/normro-ri(2)/normri)<10*eps && abs(ro(3)/normro-ri(3)/normri)<10*eps))
                    lambda=((O1(comp_idx_2)-I1(comp_idx_2)+I1(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1)-O1(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1))/(-ri(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1)+ri(comp_idx_2)));
                    xi=(I1(comp_idx_1)-O1(comp_idx_1)+lambda*ri(comp_idx_1))/ro(comp_idx_1);
                    %U1=I1+lambda*ri;
                    U1=O1+xi*ro; %xi is the distance between O1 and U1 on vector O1 to O2 normed to the total distance between O1 and O2
                    % check if point interpolated ok
                    %errU1=norm(U1-O1-lambda*ro);
                   % --Correction by Simon
                   % U1 lambda and xi better determined as intersection between plane normal to beam element and 
                   % panel. did not work for high dihedral because of strange comp_idx
                    lambda=dot((O2-I1),cut_norm)/dot((I2-I1),cut_norm); %distance of U1 and I1 normed on the vector (I2-I1)
                    U1=I1+(I2-I1)*lambda; %Intersection of Panel with current Blockface
                    xi=(U1(1)-O1(1))/ro(1);
                    % ------------------
                    % check if point is in between P1 and P2
                    %lambda=(U1(1)-I1(1))/(I2(1)-I1(1));
                    %errU1_online=norm(U1-I1-lambda*(I2-I1));
                    if lambda<1.0-100*eps && lambda>100*eps && xi<1.0+100*eps && xi>-100*eps
                        if debug==1
                            figure;
                            hold on
                            h=fill3(master_quad(1,:),master_quad(2,:),master_quad(3,:),'b');
                            set(h,'facealpha',.2);
                            h=fill3(panel(1,:),panel(2,:),panel(3,:),'r');
                            set(h,'facealpha',.2);
                            scatter3(O3(1),O3(2),O3(3),'b');
                            text(O3(1),O3(2),O3(3),'O3');
                            scatter3(O2(1),O2(2),O2(3),'b');
                            text(O2(1),O2(2),O2(3),'O2');
                            scatter3(O1(1),O1(2),O1(3),'b');
                            text(O1(1),O1(2),O1(3),'O1');
                            scatter3(I2(1),I2(2),I2(3),'r', 'filled');
                            text(I2(1),I2(2),I2(3),'I2');
                            scatter3(I1(1),I1(2),I1(3),'r', 'filled');
                            text(I1(1),I1(2),I1(3),'I1');
                            scatter3(U1(1),U1(2),U1(3),'gd');
                            text(U1(1),U1(2),U1(3),'U1');
                            grid on;
                            axis equal;
                            
                        end
                        k=k+1;
                        newpoints(:,k)=U1;
                        point_face(k)=j;
                        point_dist(k)=lambda;
%                         if errU1>1E-10
%                             %disp('geometry error');
%                         end
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
                    lambda=((O1(comp_idx_2)-I1(comp_idx_2)+I1(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1)-O1(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1))/(-ri(comp_idx_1)*ro(comp_idx_2)/ro(comp_idx_1)+ri(comp_idx_2)));
                    xi=(I1(comp_idx_1)-O1(comp_idx_1)+lambda*ri(comp_idx_1))/ro(comp_idx_1);
                    %U1=I1+lambda*ri;
                    U1=O1+xi*ro; %xi is the distance between O1 and U1 on vector O1 to O2 normed to the total distance between O1 and O2
                    % check if point interpolated ok
                  %  errU1norm(U1-O1-lambda*ro);
                   % --Correction by Simon
                   % U1 lambda and xi better determined as intersection between plane normal to beam element and 
                   % panel. did not work for high dihedral because of strange comp_idx
                    seg_norm=cross((O2-O1),(O3-O1));
                    cut_norm=cross(seg_norm,(O2-O1)); % normal vector of cutting plane
                    lambda=dot((O2-I1),cut_norm)/dot((I2-I1),cut_norm);%distance of U1 and I1 normed on the vector (I2-I1)
                    U1=I1+(I2-I1)*lambda;%Intersection of Panel with current Blockface
                    xi=(U1(1)-O1(1))/ro(1);
                    % ------------------
                    % check if point is in between P1 and P2
                    %lambda=(U1(1)-I1(1))/(I2(1)-I1(1));
                    %errU1_online=norm(U1-I1-lambda*(I2-I1));
                    if lambda<1.0-100*eps && lambda>100*eps && xi<1.0+100*eps && xi>-100*eps
                        k=k+1;
                        newpoints(:,k)=U1;
                        point_face(k)=j;
                        point_dist(k)=lambda;
%                         if errU1>1E-10
%                             %disp('geometry error');
%                         end
                    end
                end
           % end
        end
    end
end
[point_dist_s1,idx1]=sort(point_dist);
point_face_s1=point_face(idx1);
newpoints_s1=newpoints(:,idx1);

[point_face,idx2]=sort(point_face_s1);
point_dist=point_dist_s1(idx2);
newpoints=newpoints_s1(:,idx2);
if k>0
    new_polygon=zeros(3,k+4);
    r=0;
    for i=1:4 %for every face of panel
        new_polygon(:,i+r)=panel(:,i);
        for j=1:k %check for every new point if new point lies on current face of panel
            if point_face(j)==i
                r=r+1;
                new_polygon(:,i+r)=newpoints(:,r);
            end
        end
    end
    
    %%old sorting algorithm
%     new_polygon=zeros(3,k);
%     ll=1;
%     
%     tosort=[point_face;point_dist;1:k];
%     sorted=sortrows(tosort')';
%     for j=1:4 %for every panel face
%         if ll>=2
%            if ~((abs(new_polygon(1,ll-1)-panel(1,j))<100*eps) && (abs(new_polygon(2,ll-1)-panel(2,j))<100*eps) &&(abs(new_polygon(3,ll-1)-panel(3,j))<100*eps))
%                new_polygon(:,ll)=panel(:,j);
%                ll=ll+1;
%            end 
%         else
%             new_polygon(:,ll)=panel(:,j);
%             ll=ll+1;
%         end
%         for k=1:length(point_face)
%             if sorted(1,k)==j
%                if ~((abs(new_polygon(1,ll-1)-newpoints(1,sorted(3,k)))<100*eps) && (abs(new_polygon(2,ll-1)-newpoints(2,sorted(3,k)))<100*eps) &&(abs(new_polygon(3,ll-1)-newpoints(3,sorted(3,k)))<100*eps))
%                     new_polygon(:,ll)=newpoints(:,sorted(3,k));
%                     ll=ll+1;
%                 end
%             end
%         end 
%     end   
%     
else
    new_polygon=panel;
end




fact=0.0;
end



