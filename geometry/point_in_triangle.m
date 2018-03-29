
function T_in=point_in_triangle(A,B,C,P)

    
    px=P(1);
    py=P(2);
    
    p0x=A(1);
    p0y=A(2);
   
    p1x=B(1);
    p1y=B(2);
    
    p2x=C(1);
    p2y=C(2);
    
    
    Area = 1/2*(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);
    if Area==0
            px=P(1);
    py=P(3);
    
    p0x=A(1);
    p0y=A(3);
   
    p1x=B(1);
    p1y=B(3);
    
    p2x=C(1);
    p2y=C(3);
     Area = 1/2*(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);
    end
    
    
    s = 1/(2*Area)*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py);
    t = 1/(2*Area)*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py);



if (s>=-0.001) && (s<=1.001) && (t>=-0.01)  && (t<=1.001) && (s+t<=1.002)
   T_in=1;
else
    T_in=0;
end

%         u=B-A;
%     v=C-A;
%     w=P-A;
%     
%     vCrossW =cross(v,w);
%     vCrossU =cross(v,u);
%     %% TODO check limit -0.11 befor -0.05
%     if(dot(vCrossW,vCrossU) < -0.05)
%         T_in=0;
%     else
%         
%         uCrossW = cross(u, w);
%         uCrossV = cross(u, v);
% 
%         if (dot(uCrossW, uCrossV) < -0.05)
%             T_in=0;
%         else
%             denom=norm(uCrossV);
%             r = norm(vCrossW)/denom;
%             t = norm(uCrossW)/denom;
%             T_in= ((r <= 1+2E-04 )&&( t <= 1+2E-04 )&&( r + t <=1+2E-04));
%         end
%     end

end