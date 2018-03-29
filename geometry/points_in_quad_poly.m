function in=points_in_quad_poly(quad,points)
%load('debug_points_in_quad_poly.mat')
    debug=0;
    debug2=0;
    if debug==1
        close all
        figure;
        hold on
        axis equal
        h=fill3(quad(1,:),quad(2,:),quad(3,:),'g');
        set(h,'facealpha',.25)

        for i=1:size(points,2)
            scatter3(points(1,i),points(2,i),points(3,i),'b');
            text(points(1,i),points(2,i),points(3,i),['P' num2str(i)]);
        end;
    end

    red_poly=zeros(3,1);
    ll=1;
    in=zeros(size(points,2),1);
    Trias=zeros(3,3,4);
    Trias(:,:,1)=[quad(:,1),quad(:,2),quad(:,3)];
    Trias(:,:,2)=[quad(:,1),quad(:,2),quad(:,4)];
    Trias(:,:,3)=[quad(:,1),quad(:,3),quad(:,4)];
    Trias(:,:,4)=[quad(:,2),quad(:,3),quad(:,4)];
    for i=1:size(points,2)
        P=points(:,i);
        %check if point is one of the quad points
        if any((max(abs(quad-repmat(P,1,size(quad,2)))))<0.0001)
            in(i)=1;
            red_poly(:,ll)=points(:,i);
            ll=ll+1;
            if debug==1
                scatter3(P(1),P(2),P(3),'g', 'filled');
            end
        else
            inp=zeros(1,4);
            onp=zeros(1,4);
            %do for 2 triangles 1,2,3 and 1,3,4
            for k=1:size(Trias,3)

                A=Trias(:,1,k);
                B=Trias(:,2,k);
                C=Trias(:,3,k);
                %create transformation matrix which transforms the triangle to x-y
                %plane
                vx=(B-A)/norm(B-A);
                firstTriangleNormal=cross(C-A,B-A)/norm(cross(C-A,B-A));
                vz=firstTriangleNormal;
                vy=cross(vx,vz)/norm(cross(vx,vz));

                TrafoMatrix=[vx,vy,vz];
                A2=TrafoMatrix'*(A-A);
                B2=TrafoMatrix'*(B-A);
                C2=TrafoMatrix'*(C-A);
                P2=TrafoMatrix'*(P-A);
                B2=(B2)*1.01;
                C2=(C2)*1.01;
                if debug2==1
                    figure;
                    hold on
                    axis equal
                    h=fill3(quad(1,:),quad(2,:),quad(3,:),'g');
                    set(h,'facealpha',.25);
                    scatter3(A(1),A(2),A(3),'b');
                    scatter3(B(1),B(2),B(3),'b');
                    scatter3(C(1),C(2),C(3),'b');
                    scatter3(P(1),P(2),P(3),'r', 'filled');
                    text(A(1),A(2),A(3),'A');
                    text(B(1),B(2),B(3),'B');
                    text(C(1),C(2),C(3),'D');
                    text(P(1),P(2),P(3),'P');
                    line1=[A A+vx.*0.1];
                    line2=[A A+vy.*0.1];
                    line3=[A A+vz.*0.1];
                    line(line1(1,:),line1(2,:),line1(3,:));
                    line(line2(1,:),line2(2,:),line2(3,:));
                    line(line3(1,:),line3(2,:),line3(3,:));
                    text(line1(1,2),line1(2,2),line1(3,2),'vX');
                    text(line2(1,2),line2(2,2),line2(3,2),'vY');
                    text(line3(1,2),line3(2,2),line3(3,2),'vZ'); 
                    figure
                    hold on
                    axis equal
                    title('projected view normal to block')
                    scatter3(A2(1),A2(2),A2(3),'b');
                    scatter3(B2(1),B2(2),B2(3),'b');
                    scatter3(C2(1),C2(2),C2(3),'b');
                    scatter3(P2(1),P2(2),P2(3),'b+');
                    line([A2(1); B2(1)],[A2(2); B2(2)],[A2(3); B2(3)])
                    line([C2(1); B2(1)],[C2(2); B2(2)],[C2(3); B2(3)])
                    line([C2(1); A2(1)],[C2(2); A2(2)],[C2(3); A2(3)])
                    text(A2(1),A2(2),A2(3),'A');
                    text(B2(1),B2(2),B2(3),'B');
                    text(C2(1),C2(2),C2(3),'D');
                    text(P2(1),P2(2),P2(3),'P');
                end


                %check if current point lies within triangle
                [inp(k),onp(k)]=inpolygon(P2(1),P2(2),[A2(1) B2(1) C2(1)],[A2(2) B2(2) C2(2)]);
                if (inp(k)==0 && onp(k)==0)
                    [inp(k),onp(k)]=inpolygon(round(P2(1),10),round(P2(2),10),[A2(1) B2(1) C2(1)],[A2(2) B2(2) C2(2)]);
                end

            end
            if (max(inp)==0 && max(onp)==1)
                disp('bla')
            end
            if (max(inp)==1 || max(onp)==1)
                if debug==1
                    scatter3(P(1),P(2),P(3),'g', 'filled');
                end
                in(i)=1;
                red_poly(:,ll)=points(:,i);
                ll=ll+1;
            else
                if debug==1
                    scatter3(P(1),P(2),P(3),'r', 'filled');
                end
            end
        end
    end
%end