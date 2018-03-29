function in=points_in_triangle_poly(triangle,points)
red_poly=zeros(3,1);
ll=1;
for i=1:size(points,2)
    
    A=triangle(:,1);
    B=triangle(:,2);
    C=triangle(:,3);
    
    %fill3(triangle(1,:),triangle(2,:),triangle(3,:),'y')
    %plot3(points(1,i),points(2,i),points(3,i),'mx');
    
    triangle=[A B C];
     
    T1_in=point_in_triangle(A,B,C,points(:,i));

    if(T1_in==1)
        in(i)=1;
        red_poly(:,ll)=points(:,i);
        ll=ll+1;
    else
        in(i)=0;
    end
    
end
end