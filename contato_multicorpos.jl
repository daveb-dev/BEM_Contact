include("dad_contato.jl")
include("define_subregioes.jl")
include("calcula_noscontato.jl")
include("format_dad.jl")
include("calcula_centro.jl")
include("calc_fforma.jl")
include("calc_dfforma.jl")
include("calc_solfund.jl")
include("Gauss_Legendre.jl")
include("calc_ghnsing.jl")
include("calc_gsing.jl")
include("monta_GeH.jl")
include("aplica_CDC.jl")
include("gera_MatrizesGlobais.jl")
include("calc_dist.jl")
include("newton.jl")
include("verifica_contato.jl")
include("aplica_contato.jl")
include("trac_contato.jl")
include("reordena.jl")
#include("foo.jl")
using LinearAlgebra

#function contato() # For debugging

PONTOS,SEGMENTOS,MALHA,supcontato,CCSeg,E,ni,mi,cargav,w,R,tipo_problema,NOS_RES=dad_contato()
SUBREGIAO,interfaces=define_SubRegioes(SEGMENTOS,MALHA,supcontato)
noscontato,elemcontato=calcula_noscontato(interfaces,MALHA)

nnoscontato=size(noscontato,2)

NOS,ELEM,CDC,E,ni=format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg,E,ni,tipo_problema)

somalisttrac=zeros(6*size(ELEM,1),1)
somatrac=zeros(size(ELEM,1),7)
somatrac[:,1]=(1:size(ELEM,1))'
nnos=size(NOS,1)
desl=zeros(nnos,3)
nlinhas=2*nnos # Number of lines of matrix A
x0=ones(nlinhas+4*size(noscontato,2),1)
nreg=size(SUBREGIAO,1) # Number of lines in subregion i
nseg=zeros(1,size(SUBREGIAO,1)) # Number of lines in subregion i
nsegmax=size(SUBREGIAO,2) # Number of columns in matrix SUBREGIAO
nelem=zeros(1,nreg) # Ok until here

for j=1:nreg # Loop over all subregions
    i=0
    fim=1
    while (i<nsegmax && fim==1)
        seg=SUBREGIAO[j,i+1]
        if(seg==0)
            fim=0
            nseg[j]=i
        else
            nelem[j]=nelem[j]+MALHA[seg,2]
            i=i+1
        end
    end
    if(fim==1)
        nseg[j]=i
    end
end

nlinhasAA=zeros(1,nreg) # Number of lines of global A matrix with extra lines
nlinhasA=zeros(1,nreg) # Number of lines of global A matrix without extra lines
NOSt=zeros(size(NOS))
aant=0
bant=0
nelem=zeros(Int64,1,nreg)

# Creates the linear part of matrix A

ELEMloc=[]
CDCloc=[]
NOSloc=[]
G=[]
H=[]
#nelem=zeros(Int64,1,2)



A2=zeros(2068,2068)
b2=zeros(2068,1)
aant=0
bant=0

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

## Solve for multi step ##########################################################################################

npassos=1
for i=1:5*npassos
	global x0
	global somalisttrac
    println(i)
    if(i>=npassos+1)
        inilinha3=MALHA[1,2]+MALHA[2,2]+1
        fimlinha3=inilinha3+MALHA[3,2]-1
        NOS_RES=[]
        if(i==npassos+1)
            for j=inilinha3:fimlinha3
                CDC[j,2:end]=[1 cargah/npassos 0 0 1 cargah/npassos 0 0 1 cargah/npassos 0 0]
            end
#            mostra_geo[SEGMENTOS,PONTOS,ELEM,NOS] # mostra a geometria do problema
#            mostra_cdc[SEGMENTOS,PONTOS,NOS,CDC,ELEM,CCSeg,NOS_RES]
            
        elseif[i==2*npassos+1]
            for j=inilinha3:fimlinha3
                CDC[j,2:end]=[1 -cargah2/npassos 0 0 1 -cargah2/npassos 0 0 1 -cargah2/npassos 0 0]
            end
#            mostra_geo[SEGMENTOS,PONTOS,ELEM,NOS] # mostra a geometria do problema
#            mostra_cdc[SEGMENTOS,PONTOS,NOS,CDC,ELEM,CCSeg,NOS_RES]

        elseif[i==3*npassos+1]
            for j=inilinha3:fimlinha3
                CDC[j,2:end]=[1 -cargah/npassos 0 0 1 -cargah/npassos 0 0 1 -cargah/npassos 0 0]
            end
#            mostra_geo[SEGMENTOS,PONTOS,ELEM,NOS] # mostra a geometria do problema
#            mostra_cdc[SEGMENTOS,PONTOS,NOS,CDC,ELEM,CCSeg,NOS_RES]

        elseif[i==4*npassos+1]
            for j=inilinha3:fimlinha3
                CDC[j,2:end]=[1 cargah2/npassos 0 0 1 cargah2/npassos 0 0 1 cargah2/npassos 0 0]
            end
#            mostra_geo[SEGMENTOS,PONTOS,ELEM,NOS] # mostra a geometria do problema
#            mostra_cdc[SEGMENTOS,PONTOS,NOS,CDC,ELEM,CCSeg,NOS_RES]
            
        end
        aant=0
        bant=0
        for icorpo=1:2
            if(icorpo==1)
                H=H1
                G=G1
            else
                H=H2
                G=G2
            end
            seg=SUBREGIAO[icorpo,:]
            if(icorpo>1)
                nelemant=sum(nelem[1:icorpo-1])
            else
                nelemant=0
            end
            ELEMloc=ELEM[nelemant+1:nelemant+nelem[icorpo],:]
            CDCloc=CDC[nelemant+1:nelemant+nelem[icorpo],:]
            NOSloc=NOS[2*nelemant+1:2*nelemant+2*nelem[icorpo],:]
            A,B,DESL_PRloc,b=aplica_CDC(ELEMloc,NOSloc,H,G,CDCloc,NOS_RES)
            AA,bb=gera_MatrizesGlobais(A,B,b,ELEMloc,nnoscontato)
            a,b = size(AA)
            ndfo = aant + 1
            ndff = a + aant
            ndco = bant + 1
            ndcf = b + bant
            aant = ndff
            bant = ndcf
            A2[ndfo:ndff,ndco:ndcf] = AA # Joga as matrizes locais na global
            b2[ndfo:ndff,1] = bb
            if(icorpo==1)
                DESL_PR1=DESL_PRloc
            else
                DESL_PR2=DESL_PRloc
            end
        end
        milocal=mi
    else
        milocal=mi
    end
    NOSt[:,2]=NOS[:,2]+desl[:,2] # Posição dos nós + desl
    NOSt[:,3]=NOS[:,3]+desl[:,3] # Posição dos nós + desl
    # Cálculo da distância h entre os nós na região de provável
    # contato
    h=calc_dist(NOSt,noscontato)
    println("Newton")
    erromax = 2e-10
    x,contato=newton(x0,A2,b2,h,milocal,somalisttrac,nnoscontato,nlinhasA,nlinhasAA,noscontato,erromax)
    if(i==2*npassos||i==npassos||i==4*npassos) 
        x0=ones(nlinhas+4*size(noscontato,2),1)
    else
        x0=x
    end
    for icorpo=1:2
        if(icorpo==1)
            ilinha=0
            ielem=0
            ino=0
            DESL_PR=DESL_PR1
        else
            ilinha=nlinhasAA[icorpo-1]
            ino=nelem[icorpo-1]*2
            ielem=nelem[icorpo-1]
           DESL_PR=DESL_PR2
        end
listdesllocal,listtra,desllocal,traclocal = reordena[x[1+ilinha:nlinhasAA[icorpo]+ilinha], CDC[1+ielem:nelem[icorpo]+ielem,:],DESL_PR] 

somalisttrac[6*ielem+1:6*(ielem+nelem[icorpo])]= somalisttrac[6*ielem+1:6*(ielem+nelem[icorpo])]+listtra
       
somatrac[1+ielem:nelem[icorpo]+ielem,:]=somatrac[1+ielem:nelem[icorpo]+ielem,:]+traclocal
        desl[1+ino:ino+nelem[icorpo]*2,2:end]=desl[1+ino:ino+nelem[icorpo]*2,2:end]+desllocal[:,2:end]
        listdesl[2*(1+ino)-1:2*(ino+nelem[icorpo]*2),:]=listdesl[2*(1+ino)-1:2*(ino+nelem[icorpo]*2),:]+listdesllocal
    end
    nlisttra=size(somalisttrac,1)
    listtra,trac=trac_contato(contato,x,nlisttra,noscontato,nlinhasA,nlinhasAA,nelem)
    somalisttrac=somalisttrac+listtra
    somatrac[:,2:7]=somatrac[:,2:7]+trac[:,2:7]    
end

#end # For debugging
