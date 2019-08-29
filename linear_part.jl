# Liner part of matrix

function linear_part()

A2=zeros(2068,2068)
b2=zeros(2068,1)

for icorpo=1:2
	global aant
	global bant	
	
    if(icorpo==2)
    global NOS_RES=[]
    #NOS_RES=[]
    end
    seg=SUBREGIAO[icorpo,:]
    nelem[icorpo]=sum(MALHA[seg,2])
    if(icorpo>1)
        nelemant=sum(nelem[1:icorpo-1])
    else
        nelemant=0
    end
    global ELEMloc=ELEM[nelemant+1:nelemant+nelem[icorpo],:]
    global CDCloc=CDC[nelemant+1:nelemant+nelem[icorpo],:]
    global NOSloc=NOS[2*nelemant+1:2*nelemant+2*nelem[icorpo],:]
    G,H=monta_GeH(ELEMloc,NOSloc,E,ni)
    #ELEMloc=ELEM[nelemant+1:nelemant+nelem[icorpo],:]
    #CDCloc=CDC[nelemant+1:nelemant+nelem[icorpo],:]
    #NOSloc=NOS[2*nelemant+1:2*nelemant+2*nelem[icorpo],:]
    #G,H=monta_GeH(ELEMloc,NOSloc,E,ni) # Assembles H and G matrices
    # Aplica as condições de contorno considerando a região de slip como se
    # fosse região de contato sem atrito; ou seja; tx é conhecido e igual a
    # zero.
    if(icorpo==1)
       H1=H
       G1=G
    else
       H2=H
       G2=G
    end

    A,B,DESL_PRloc,b=aplica_CDC(ELEMloc,NOSloc,H,G,CDCloc,NOS_RES)
    AA,bb=gera_MatrizesGlobais(A,B,b,ELEMloc,nnoscontato)
    a,b = size(AA)
    ndfo = aant + 1
    ndff = a + aant
    ndco = bant + 1
    ndcf = b + bant
    aant = ndff
    bant = ndcf
    #global A2[ndfo:ndff,ndco:ndcf] = AA # Joga as matrizes locais na 
    #global b2[ndfo:ndff,1] = bb
    A2[ndfo:ndff,ndco:ndcf] = AA
    b2[ndfo:ndff,1] = bb
    if(icorpo==1)
        DESL_PR1=DESL_PRloc
    else
        DESL_PR2=DESL_PRloc
    end
    nlinhasAA[icorpo]=size(AA,1) # Número de linhas da matriz AA
    nlinhasA[icorpo]=size(A,1) # Número de linhas da matriz A	
end

sum_nlinhasA=sum(nlinhasA)
sum_nlinhasA=floor(Int,sum_nlinhasA)
listdesl=zeros(sum_nlinhasA,1)

return A2,b2
end
