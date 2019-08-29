function calc_ghnsing(x1,y1,x2,y2,x3,y3,x0,y0,E,ni,qsi,w)

# Function que calcula as matrizes g e h locais quando o ponto fonte n�o
# pertence ao elemento [integra��o regular]

npg=length(qsi) # N�mero de pontos de Gauss

h=zeros(2,6) # Inicializa a matrizcalc h do elemento
g=zeros(2,6) # Inicializa a matriz g do elemento

for i=1:npg
    N1,N2,N3=calc_fforma(qsi[i]) # Calcula as fun��es de forma
    N=[N1 0 N2 0 N3 0;0 N1 0 N2 0 N3]

    x=N1*x1+N2*x2+N3*x3 # Calcula a coordenada x do ponto de integra��o
    y=N1*y1+N2*y2+N3*y3 # Calcula a coordenada y do ponto de integra��o
    
    dN1dqsi,dN2dqsi,dN3dqsi=calc_dfforma(qsi[i]) # Calcula as 
                                           # derivadas das fun��es de forma
    
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2)
    
    sx=1/dgamadqsi*dxdqsi # Component x do vetor tangente
    sy=1/dgamadqsi*dydqsi # Componente y do vetor tangente
    
    nx=sy; # Componente x do vetor normal
    ny=-sx; # Componente y do vetor normal
    
    uast,tast=calc_solfund(x,y,x0,y0,nx,ny,E,ni) # Calcula as solu��es fundamentais
    h=h+tast*N*dgamadqsi*w[i] # Integra��o da matriz h see error
    g=g+uast*N*dgamadqsi*w[i] # Integra��o da matriz g
end
return g,h
end
