function calc_gsing(x1,y1,x2,y2,x3,y3,qsi0,E,ni)

# Pontos de Gauss Padr�o
qsi = [-0.97390652851717 0.97390652851717 -0.86506336668898 0.86506336668898 -0.67940956829902 0.67940956829902 -0.43339539412925 0.43339539412925 -0.14887433898163 0.14887433898163]
# Pesos de Gauss Padr�o
w = [ 0.06667134430868 0.06667134430868 0.14945134915058 0.14945134915058 0.21908636251598 0.21908636251598 0.26926671931000 0.26926671931000 0.29552422471475 0.29552422471475]
# Pontos de Gauss logar�tmico
eta = [0.90425944e-2 0.53971054e-1 0.13531134 0.24705169 0.38021171 0.52379159 0.66577472 0.79419019 0.89816102 0.96884798]

# Pesos de Gauss Logar�tmico
rho = [0.12095474 0.18636310 0.19566066 0.17357723 0.13569597 0.93647084e-1 0.55787938e-1 0.27159893e-1 0.95151992e-2 0.16381586e-2]

npgr=length(qsi) # N�mero de pontos de Gauss padr�o [regular]
npgs=length(eta) # N�mero de pontos de Gauss logar�tmico [singular]

gs=zeros(2,2) # Inicializa a integral singular
gns=zeros(2,2) # Inicializa a integral n�o singular

GEL=E÷(2*(1+ni)) # M�dulo de elasticidade transversal
for i=1:npgr # La�o sobre os pontos de integra��o regular
    N1,N2,N3=calc_fforma(qsi[i]) # Fun��es de forma
    if(qsi0==-1)
        deno=1
        x0=x1
        y0=y1
        ff = [N1 0; # Vetor de fun��es de forma
            0 N1]
    else()
        deno=2
        x0=x1
        y0=y1
        ff = [N2 0; # Vetor de fun��es de forma
            0 N2]
        
    end
    
    dN1dqsi,dN2dqsi,dN3dqsi=calc_dfforma(qsi[i]) # Derivadas das fun��es de forma
    
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2)
    
    x=N1*x1+N2*x2+N3*x3; # Calcula a coordenada x do ponto de integra��o
    y=N1*y1+N2*y2+N3*y3; # Calcula a coordenada y do ponto de integra��o
    r=sqrt((x-x0)^2+(y-y0)^2)
    rx=(x-x0)/r
    ry=(y-y0)/r
    rd=[rx^2 rx*ry; rx*ry ry^2]
    l1 = ((qsi[i]+qsi0)*(x1-2*x2+x3)+x3-x1)/deno
    l2 = ((qsi[i]+qsi0)*(y1-2*y2+y3)+y3-y1)/deno
    
    lnrns = -log(sqrt(l1^2+l2^2))
    
    # Parte n�o singular da solu��o fundamental de deslocamento [u_est_ns[2x2]]
    
    u_est_ns = ((3-4*ni)*[lnrns 0 ; 0 lnrns] + rd)/(8*pi*GEL*(1-ni))
    
    # Itera��o de soma da integra��o n�o singular
    gns = gns + u_est_ns*ff*dgamadqsi*w[i]
    
end

# Integra��o singular
for i=1:npgs
    if(qsi0==-1) # O ponto fonte � o primeiro (ou o terceiro, uma vez que
        # a integral � a mesma)
        qsil=2*eta[i]-1
        N1,N2,N3=calc_fforma(qsil) # Fun��es de forma
        dqsideta=2
        ff=[N1 0;0 N1]
    else()
        qsil=eta[i]; # O ponto fonte � o do centro
        N1,N2,N3=calc_fforma(qsil) # Fun��es de forma
        dqsideta=1
        ff=[N2 0;0 N2]
    end
    dN1dqsi,dN2dqsi,dN3dqsi=calc_dfforma(qsil) #Derivada das fun��es de forma
    dxdqsi=dN1dqsi*x1+dN2dqsi*x2+dN3dqsi*x3
    dydqsi=dN1dqsi*y1+dN2dqsi*y2+dN3dqsi*y3
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2)
    if(qsi0==-1)
        gs=gs+(3-4*ni)/(8*pi*GEL*(1-ni))*ff*dgamadqsi*rho[i]*dqsideta
    else()
        gs=gs+2*(3-4*ni)/(8*pi*GEL*(1-ni))*ff*dgamadqsi*rho[i]*dqsideta
    end
end
g=gs+gns
return g
end
