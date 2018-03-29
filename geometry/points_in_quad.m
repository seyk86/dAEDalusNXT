function in=points_in_quad(quad,points)

for i=1:size(points,2)
    
    A=quad(:,1);
    B=quad(:,2);
    C=quad(:,4);
    
    triangle=[A B C];
    
%     fill3(triangle(1,:),triangle(2,:),triangle(3,:),'y')
%     plot3(points(1,i),points(2,i),points(3,i),'bx');
   
    T1_in=point_in_triangle(A,B,C,points(:,i));
    
    A=quad(:,2);
    B=quad(:,3);
    C=quad(:,4);
    triangle=[A B C];
%     fill3(triangle(1,:),triangle(2,:),triangle(3,:),'y')
%     plot3(points(1,i),points(2,i),points(3,i),'bx');
%     
    T2_in=point_in_triangle(A,B,C,points(:,i));
    
    if(T1_in==1 || T2_in==1)
        in(i)=1;
    else
        in(i)=0;
    end
    
end
end
