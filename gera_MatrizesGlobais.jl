function gera_MatrizesGlobais(A,B,b,ELEMloc,nnoscontato)

nlinhas=size(A,2) # N�mero de linhas da matriz A
A2=zeros(nlinhas+2*nnoscontato,nlinhas+2*nnoscontato) # Inicializa��o da
# matriz A2 que � composta pela matriz A mais as colunas da matriz B
# referentes aos n�s da regi�o de contato e mais as linhas referentes �s
# linhas das condi��es de contato
A2[1:nlinhas,1:nlinhas]=A # Insere a matriz A dentro de A2
b2=[b;zeros(2*nnoscontato,1)] # Inicializa��o do vetor b2 que � formado
# pelo vetor b mais algumas linhas referentes �s condi��es de contato
primeiro=ELEMloc[1,2]
no_contato::Int32=0

## Constru��o das matrizes A2 e B sem as condi��es de contato

for k=1:(nnoscontato-1)÷2 # Percorre todos os elementos que podem entrar
    #     #    em contato
    for no=1:3 # Percorre os 3 n�s do elemento el
        no_contato=ELEMloc[k,no+1]-primeiro+1 # n�mero global do n�
        postx=nlinhas+2*no_contato-1  # Posi��o da coluna da matriz A2
        # referente � for�a de superf�cie na dire��o x
        posty=nlinhas+2*no_contato  # Posi��o da coluna da matriz A2
        # referente � for�a de superf�cie na dire��o y
        poslinhasA=[postx posty] # Colunas de A que receber�o colunas
        # de B
        poscolunasB=6*(k-1)+2*(no-1)+1:6*(k-1)+2*(no-1)+2
        colunasB=B[:,poscolunasB] # colunas de B que ir�o para A2

        A2[1:nlinhas,postx:posty] = A2[1:nlinhas,postx:posty]-colunasB
        B[:,poscolunasB]=zeros(nlinhas,2) # Zera as colunas de B
    end
end
return A2,b2
end
