function calcula_centro(x1,y1,x2,y2,raio)
    # Compute the center of an arc given two points and the radius

    xm=(x1+x2)/2
    ym=(y1+y2)/2
    b=√((x2-x1)^2+(y2-y1)^2)
    t1=(x2-x1)/b
    t2=(y2-y1)/b
    n1=t2
    n2=-t1
    h=√(abs(raio^2-(b/2)^2))
    if(raio>0)
        if(n1==0)
            xc=xm
            yc=ym-n2/abs(n2)*h
        else
            xc=-n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
            yc=n2/n1*(xc-xm)+ym
        end
    else
        if(n1==0)
            xc=xm
            yc=ym+n2/abs(n2)*h
        else
            xc=n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
            yc=n2/n1*(xc-xm)+ym
        end
    end
    xc,yc
end

function calcula_arco(x1,y1,x2,y2,xc,yc,raio)
# Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
# horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
# and the horizontal direction

	dx1 = x1 - xc; dy1 = y1 - yc
	dx2 = x2 - xc; dy2 = y2 - yc
	tet1=atan(dy1,dx1);
	a=[dx1,dy1,0];
	b=[dx2,dy2,0];
	angle = atan(norm(cross(a,b)),dot(a,b));
	if(raio>0)
	    tet2=tet1+angle;
	else
	    tet2=tet1-angle;
	end
	tet1,tet2
end
