function [ area ] = planar_polygon_area(polygon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%create normal to plane out of first 3 points
v12=polygon(:,2)-polygon(:,1);
v13=polygon(:,3)-polygon(:,1);
vn=cross(v12,v13);
if norm(vn)==0
    if size(polygon,2)>3
        v14=polygon(:,4)-polygon(:,1);
        vn=cross(v12,v14);
        %transform polygon to x-y plane
        vx=v12/norm(v12);
        vz=vn;
        vy=cross(vx,vn);
        TrafoMatrix=[vx,vy,vz];

        transformedpolygon=TrafoMatrix'*polygon;

        %use polyarea
        area=polyarea(transformedpolygon(1,:),transformedpolygon(2,:));
    else
        area=0;
    end
else
    vn=vn/norm(vn);

    %check if all other points lie within this plane
    % if size(polygon,2)>3
    %     for i=1:size(polygon,2)
    %         v1i=polygon(:,i)-polygon(:,1);
    %         if dot(v1i,vn)>100*eps
    %            disp('warning, points do not lie within one plane');
    %         end
    %     end
    % end
    %transform polygon to x-y plane
    vx=v12/norm(v12);
    vz=vn;
    vy=cross(vx,vn);
    TrafoMatrix=[vx,vy,vz];

    transformedpolygon=TrafoMatrix'*polygon;

    %use polyarea
    area=polyarea(transformedpolygon(1,:),transformedpolygon(2,:));
end
end

%OLD WAY
% 
% area=0;
% coord=3;
% 
% % TODO for 3D
% 
% if abs(polygon(2,1))<100*eps
%     idx1=1;
%     idx2=3;
%     nidx=2;
% else
%     idx1=1;
%     idx2=2;
%     nidx=3;
%     
% end
% 
% 
% n=cross(polygon(:,1)-polygon(:,end-1),polygon(:,2)-polygon(:,3));
% if (sum(n)<=1E-10) && (sum(n)>=-1E-10)
%     area=0;
% else
%     n=n/norm(n);
%     area=polyarea(polygon(idx1,:),polygon(idx2,:));
%     if isinf(area)
%         area=0;
%     else
%         if abs(n(nidx))>1e-15 
%             area=area/n(nidx);
%         else
%             area=0;
%         end
%     end
% end
% %end
% 
% 