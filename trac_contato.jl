function trac_contato(contato,x,nlisttra,noscontato,nlinhasA,nlinhasAA,nelem)

## For�as de superf�cie referentes � regi�o de poss�vel contato
nnoscontato=size(contato,1)
trac=zeros(nlisttra,1)
ndof=nlinhasA[1]
nelemtotal=sum(nelem)
T=zeros(nelemtotal,7)
ndof2=nlinhasA[2]+nlinhasAA[1]
for k=1:nnoscontato # Percorre todos os elementos que podem entrar
    if(k==50)
#        disp(' ')
    end
    noglobal=k # n�mero do n� em contato
    noglobal2=noscontato[2,k]-noscontato[2,nnoscontato]+1 # n�mero do n� em contato
    el=ceil(k/2)
    el2=floor(noglobal2/2)
    no=k-2*el+2
    no2=noglobal2-2*el2+2
    posx=ndof+2*noglobal-1
    posy=ndof+2*noglobal
    posx2=ndof2+2*noglobal2-1
    posy2=ndof2+2*noglobal2
    valorx=x[posx]
    valory=x[posy]
    valorx2=x[posx2]
    valory2=x[posy2]
    T[el,2*(no-1)+2]=valorx
    T[el,2*(no-1)+3]=valory
    indx=(6*el+2*no-8)+1
    indy=(6*el+2*no-8)+2
    trac[indx]=valorx
    trac[indy]=valory
    if(el2!=0)
        T[el2+nelem[1],2*(no2-1)+2]=valorx2
        T[el2+nelem[1],2*(no2-1)+3]=valory2
        indx2=(6*(el2+nelem[1])+2*no2-8)+1
        indy2=(6*(el2+nelem[1])+2*no2-8)+2
        trac[indx2]=valorx2
        trac[indy2]=valory2
    else
        T[nelem[1]+1,2]=valorx2
        T[nelem[1]+1,3]=valory2
        indx2=(6*(1+nelem[1])+2*1-8)+1
        indy2=(6*(1+nelem[1])+2*1-8)+2
        trac[indx2]=valorx2
        trac[indy2]=valory2        
    end
    if(no==1&&el!=1) # repete os valores de n�s compartilhados
        T[el-1,6]=valorx # tx do n� 3 do elemento el-1
        T[el-1,7]=valory # ty do n� 3 do elemento el-1
        T[nelem[1]+el2+1,2]=valorx2
        T[nelem[1]+el2+1,3]=valory2
        indx=(6*(el-1)+2*3-8)+1
        indy=(6*(el-1)+2*3-8)+2
        indx2=(6*(el2+nelem[1]+1)+2*1-8)+1
        indy2=(6*(el2+nelem[1]+1)+2*1-8)+2
        trac[indx]=valorx
        trac[indy]=valory
        trac[indx2]=valorx2
        trac[indy2]=valory2
    end
end
return trac,T
end
