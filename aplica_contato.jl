function aplica_contato(A2,b2,h,mi,somalisttrac,contato,nlinhasA,nlinhasAA,noscontato)
nnoscontato=size(contato,1) #
# Cria a parte não linear da matriz A2
## Imposição das condições de contato na matriz A2 e construção do vetor b2
for k=1:nnoscontato # Percorre todos os elementos que podem entrar
    tipocontato=contato[k,2] # tipo da condição de contato
    no_contato1=noscontato[1,k] # número do nó em contato
    no_contato2=noscontato[2,k]-noscontato[2,nnoscontato]+1 # número do nó em contato
    posux1::Int32=2*no_contato1-1 # Posição da coluna da matriz A2 referente
    # ao deslocamento na direção x
    posux2::Int32=2*no_contato2-1+nlinhasAA[1] # Posição da coluna da matriz A2 referente
    # ao deslocamento na direção x
    posuy1::Int32=2*no_contato1 # Posição da coluna da matriz A2 referente
    # ao deslocamento na direção y
    posuy2::Int32=2*no_contato2+nlinhasAA[1] # Posição da coluna da matriz A2 referente
    # ao deslocamento na direção y
    postx1::Int32=nlinhasA[1]+2*no_contato1-1  # Posição da coluna da matriz A2
    # referente à força de superficie na direção x
    postx2::Int32=nlinhasA[2]+nlinhasAA[1]+2*no_contato2-1  # Posição da coluna da matriz A2
    # referente à forçaa de superfície na direção x
    posty1::Int32=nlinhasA[1]+2*no_contato1  # Posição da coluna da matriz A2
    posty2::Int32=nlinhasAA[1]+nlinhasA[2]+2*no_contato2  # Posição da coluna da matriz A2
    el=ceil(no_contato1/2)
    no=no_contato1-2*el+2
#    poscolunasB=6*(el-1)+2*(no-1)+1:6*(el-1)+2*(no-1)+2
    poscolunasB=Int32[6*(el-1)+2*(no-1)+1 6*(el-1)+2*(no-1)+2]
    if(tipocontato==1) # Zona livre de contato
        A2[postx1,postx1] = 1; # tx1=0
        A2[posty1,posty1] = 1; # ty1=0
        A2[postx2,postx2] = 1; # tx2=0
        A2[posty2,posty2] = 1; # ty2=0
    end
    if((tipocontato==2||tipocontato==-2)) # Zona em deslizamento
        sinal=sign(tipocontato)
        A2[postx1,postx1] = 1 # tx1=-tx2
        A2[postx1,postx2] = 1 # tx1=-tx2
        A2[posty1,postx1] = 1 # tx1 = -sinal[ux]*mi abs(ty1)
        A2[posty1,posty1] = sinal*mi # tx1 = -sinal[ux]*mi abs(ty1)
        A2[postx2,posuy1] = 1 # uy1-uy2=gap
        A2[postx2,posuy2] = -1 # uy1-uy2=gap
        A2[posty2,posty1] = 1 # ty1 = -ty2
        A2[posty2,posty2] = 1 # ty1 = -ty2
        b2[posty1]=b2[posty1]-somalisttrac[poscolunasB[2]]*sinal*mi-somalisttrac[poscolunasB[1]] # tx1= mi ty1 # error is here
        b2[postx2]=b2[postx2]+h[no_contato1] # uy1-uy2=gap
    end
    if(tipocontato==3) # Zona em adesão
        A2[postx1,postx1] = 1 # tx1=-tx2
        A2[postx1,postx2] = 1 # tx1=-tx2
        A2[posty1,posux1] = 1 # ux1=ux2
        A2[posty1,posux2] = -1 # ux1=ux2
        A2[postx2,posuy1] = 1 # uy1-uy2=gap
        A2[postx2,posuy2] = -1 # uy1-uy2=gap
        A2[posty2,posty1] = 1 # ty1=-ty2
        A2[posty2,posty2] = 1 # ty1=-ty2
        b2[postx2]=b2[postx2]+h[no_contato1] # uy1-uy2=gap
    end
end
return A2,b2,contato
end
