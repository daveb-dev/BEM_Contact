function verifica_contato(x,h,mi,nlinhasA,nlinhasAA,nnoscontato,somalisttrac,noscontato)
# Verifica as condições de contato entre os nós da região de contato
contato=zeros(nnoscontato,2)

for k=1:nnoscontato # Percorre todos os elementos que podem entrar em contato
    no_contato1=noscontato[1,k]
    posux1::Int32=2*no_contato1-1; # Posição de ux do nó k no vetor x
    posuy1::Int32=2*no_contato1; # Posição de uy do nó k no vetor x
    postx1::Int32=nlinhasA[1]+2*no_contato1-1; # Posição de tx no vetor x
    posty1::Int32=nlinhasA[1]+2*no_contato1; # Posição de ty no vetor x
    no_contato2=noscontato[2,k]-noscontato[2,nnoscontato]+1; # número do nó em contato
    posuy2::Int32=2*no_contato2+nlinhasAA[1]; # Posição da coluna da matriz A2 referente
    posux2::Int32=2*no_contato2+nlinhasAA[1]-1; # Posição da coluna da matriz A2 referente
    postx2::Int32=nlinhasA[1]+nlinhasAA[1]+2*no_contato2-1;  # Posição da coluna da matriz A2
    # ao deslocamento na direção y
    el=ceil(no_contato1/2)
    no=no_contato1-2*el+2
    poscolunasB=Int32[6*(el-1)+2*(no-1)+1 6*(el-1)+2*(no-1)+2]
    tx = x[postx1]+ somalisttrac[poscolunasB[1]]
    ty = x[posty1]+ somalisttrac[poscolunasB[2]]
    ux1 = x[posux1]
    uy1 = x[posuy1]
    uy2 = x[posuy2]
    gap=h[no_contato1]
    if(abs(gap)<1e-8)
        gap=0
    end
    sinal=sign(ux1)
    if(sign(ux1)==0)
        sinal=1
    end
    if (-uy1+uy2<-gap && (abs(-uy1+uy2+gap)>abs(gap)*1e-8 && gap!=0)) ||ty<0
        contato[no_contato1,1]=no_contato1; #
        contato[no_contato1,2]=1; # Nó livre
    else
        if (abs(tx)>=mi*abs(ty)&&sign(tx/ty)==-sinal)||abs(abs(tx)-mi*abs(ty))<abs(1e-8*tx)||mi==0 # O nó desliza
            contato[no_contato1,1]=no_contato1
            contato[no_contato1,2]=sinal*2; # tipo de contato = deslizamento
            x[postx1]=(-sinal*mi*ty-somalisttrac[poscolunasB[1]])
            x[postx2]=-x[postx1]
        else # nó em adesão
            contato[no_contato1,1]=no_contato1
            contato[no_contato1,2]=3; # tipo de contato = adesão
            x[posux2]=x[posux1]
        end
    end
end
return contato
end

