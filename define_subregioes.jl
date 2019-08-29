function define_SubRegioes(SEGMENTOS,MALHA,supcontato)


# Subregions definition
#--------------------------------------------------------------------------
# REGION = [ linha 1, linha 2, linha 3, linha 4 ]
#--------------------------------------------------------------------------
nlinhas=length(SEGMENTOS[:,2])
j=1
k=1
jj=1
kk=1
SUBREGIAO=zeros(Int8,2,4)
ponto_ini=0
for i=1:nlinhas
    if(k==1)
        ponto_ini=SEGMENTOS[i,2]
    end
    p1=SEGMENTOS[i,2]
    p2=SEGMENTOS[i,3]
    SUBREGIAO[j,k]=i
    k=k+1
    if p2==ponto_ini
        j=j+1
        k=1
    end
end

interfaces=zeros(Int64,1,6)
nreg=length(SUBREGIAO[:,1]);  # Number of subregions
nlreg=length(SUBREGIAO[1,:]); # Number of lines in the subregion
nsupcontato=size(supcontato,1)
# Interfaces = [num. da interface, subreg 1, subreg 2, segmento 1, segmento 2, nnos]
for i=1:nreg # Region in analysis [fixa]
    for j=1:nsupcontato
        interfaces[j,1]=j
        for m=1:2
            isup=supcontato[j,m]
            for k=1:nlreg
                seg=SUBREGIAO[i,k]
                if(isup==seg)
                    interfaces[j,m+1]=i
                    interfaces[j,m+3]=isup
                    nelem=MALHA[isup,2]; # Discretization in the same line for two regions
                    nnos=2*nelem+1; # Constant element
                    interfaces[j,6]=nnos
                end
            end
        end
    end
end
return SUBREGIAO,interfaces
end
