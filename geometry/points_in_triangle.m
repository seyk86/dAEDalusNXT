function in=points_in_triangle(triangle,points)

for i=1:size(points,2)
    
    A=triangle(:,1);
    B=triangle(:,2);
    C=triangle(:,3);
    
    triangle=[A B C];

    T1_in=point_in_triangle(A,B,C,points(:,i));
    
    if(T1_in==1)
        in(i)=1;
    else
        in(i)=0;
    end
    
end
end
