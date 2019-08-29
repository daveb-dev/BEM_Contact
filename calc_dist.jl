function calc_dist(NOS,noscontato)
# Calcula a dist�ncia h entre os n�s na regi�o de prov�vel
# contato
ini=noscontato[1,1] # N�mero do primeiro elemento que pode entrar em contato

fim=noscontato[1,end] # N�mero do �ltimo elemento que pode entrar em contato
nnos=size(NOS,1) # N�mero total de n�s
h=zeros(nnos,1) # Cont�m a dist�ncia inicial dos n�s que podem entrar em
# contato. Para a regi�o onde n�o existe possibilidade
# de contato; h=0.
icount=0
for i=ini:fim
    icount=icount+1
    ino2=noscontato[2,icount]
    y=NOS[i,3]-NOS[ino2,3]
    h[i,1]=-y
end
return h
end
