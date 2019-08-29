function monta_GeH(ELEMloc,NOSloc,E,ni)

# Function que monta as matrizes G e H globais

n_el = length(ELEMloc[:,1])	# N�mero total de elementos
n_nos = length(NOSloc[:,1])	# N�mero total de n�s
primeiro=ELEMloc[1,2]
n_pontos_Gauss=10
qsi,w=Gauss_Legendre(-1,1,n_pontos_Gauss)
no1::Int32=0
no2::Int32=0
no3::Int32=0
noglobal::Int32=0
#println("ELEM=$ELEM")
#println("NOS=$NOS")

# Inicializa��o das matrizes H e G
H = zeros(2*n_nos,2*n_nos)
G = zeros(2*n_nos,6*n_el)

for i = 1 : n_nos	# Percorre os pontos fontes
    
    # Coordenadas [xf,yf] dos pontos fontes
    x0 = NOSloc[i,2]
    y0 = NOSloc[i,3]
    
    for j = 1 : n_el	# Percorre os elementos
        
        # Numera��o dos tr�s n�s do elemento i
        no1 = ELEMloc[j,2]-primeiro+1
        no2 = ELEMloc[j,3]-primeiro+1
        no3 = ELEMloc[j,4]-primeiro+1
        
        # Coordenadas dos tr�s n�s [x1,y1,x2,y2,x3,y3]
        x1 = NOSloc[no1,2];	y1 = NOSloc[no1,3]
        x2 = NOSloc[no2,2];	y2 = NOSloc[no2,3]
        x3 = NOSloc[no3,2];	y3 = NOSloc[no3,3]
        
        # C�lculo das submatrizes h_el e g_el
        
        # Integra��o n�o singular
        g,h = calc_ghnsing(x1,y1,x2,y2,x3,y3,x0,y0,E,ni,qsi,w)
        
        if ((i == no1) || (i == no3))  # O n� j pertence ao elemento i
            gsing=calc_gsing(x1,y1,x2,y2,x3,y3,-1,E,ni) # Integra��o 
                      #     singular com ponto fonte na extremidade 
                      #     do elemento qsi0=-1 ou qsi0=1)
            if i==no1
                g[:,1:2]=gsing
            else
                g[:,5:6]=gsing
            end
        elseif i==no2
            gsing=calc_gsing(x1,y1,x2,y2,x3,y3,0,E,ni)# Integra��o 
                     #     singular com ponto fonte no meio 
                     #     do elemento [qsi0 = 0]
            g[:,3:4]=gsing
        end
        for nolocal = 1 : 3
            noglobal = ELEMloc[j,nolocal+1]-primeiro+1 #�ndice da matriz global H
            H[2*i-1:2*i,2*noglobal-1:2*noglobal]=H[2*i-1:2*i,2*noglobal-1:2*noglobal]+h[:,2*nolocal-1:2*nolocal]
        end
        G[2*i-1:2*i,6*j-5:6*j] = g
    end
end

# Calculo dos termos da diagonal da matriz H [considera�ao de corpo a
# temperatura constante]

for m = 1 : n_nos
    H[2*m-1:2*m,2*m-1:2*m] .= 0 #zera a diagonal principal
    for n = 1 : n_nos
        if n != m
            H[2*m-1:2*m,2*m-1:2*m]=H[2*m-1:2*m,2*m-1:2*m]-H[2*m-1:2*m,2*n-1:2*n]
        end
    end
end 
return G,H
end 
