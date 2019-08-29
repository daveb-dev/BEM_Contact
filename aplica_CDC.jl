function aplica_CDC(ELEMloc,NOSloc,A,B,CDCloc,NOS_RES)
# Troca as colunas das matrizes A e B conforme as condi��es de contorno
#ELEMloc,NOSloc,H,G,CDCloc,NOS_RES

#A=H
#B=G
n_el=size(CDCloc,1)
nnos=size(NOSloc,1)
todos_valores=zeros(1,6*n_el)
nos_elem=zeros(nnos,5)
cont=0
primeiro=ELEMloc[1,2]
no1::Int32=0
no2::Int32=0
no3::Int32=0
no_global::Int32=0
ind_A::Int32=0
ind_B::Int32=0

DESL_PR=zeros(1,6)
desl_aux=zeros(1,6)
GEL=1
for el=1:n_el
    no3=ELEMloc[el,4]-primeiro+1
    for no=1:3
        no_global=ELEMloc[el,no+1]-primeiro+1
        if(no==1||no==2)
            nos_elem[no_global,1]=no_global
            nos_elem[no_global,2]=el
            nos_elem[no_global,3]=no
            if(no==1)
                nos_elem[no3,4]=el
                nos_elem[no3,5]=3
            end
        end
        for gdl=1:2
            tipoCDC=CDCloc[el,4*no-3+2*gdl-1]
            valorCDC=CDCloc[el,4*no-2+2*gdl-1]
            todos_valores[6*el+2*no+gdl-8]=valorCDC
            if(tipoCDC==0)
                compartilha=0
                i=1
                while (compartilha==false&&i<=cont&&cont>0)
                    if(no_global==DESL_PR[i,1]&&gdl==DESL_PR[i,2])
                        compartilha=1
                        DESL_PR[i,5]=el
                        DESL_PR[i,6]=no
			desl_aux[5]=el
                        desl_aux[6]=no
                    end
                    i=i+1
                end
                if(compartilha==false)
                    cont=cont+1
                    DESL_PR[cont,1]=no_global
                    DESL_PR[cont,2]=gdl
                    DESL_PR[cont,3]=el
                    DESL_PR[cont,4]=no
		    desl_aux[1]=no_global
                    desl_aux[2]=gdl
                    desl_aux[3]=el
                    desl_aux[4]=no
              	    DESL_PR=vcat(DESL_PR,desl_aux)		    	    
       		end		
            end
        end
    end
end

if(~isempty(NOS_RES))
    n_nosres=length(NOS_RES[:,1])
	n_nosres=1
    for i=1:n_nosres
        num_no=NOS_RES[i,1]
        dire=NOS_RES[i,2]
        num_elem1=nos_elem[num_no,2]
        num_loc1=nos_elem[num_no,3]
        num_elem2=nos_elem[num_no,4]
        num_loc2=nos_elem[num_no,5]
        cont=cont+1
        DESL_PR[cont,:]=[num_no,dire,num_elem1,num_loc1,num_elem2,num_loc2]
    end
end

if(DESL_PR[1,1]!=0)
    n_des_pr = size(DESL_PR[:,1],1); # Number of restricted displacements
else
    n_des_pr = 0
end

for de=1 : n_des_pr # for over the restricted nodes
    n_no=DESL_PR[de,1]
    n_gdl=DESL_PR[de,2]
    n_el=DESL_PR[de,3]
    n_no_loc=DESL_PR[de,4]
    ind_A=(2*n_no-2)+n_gdl
    ind_B=(6*n_el+2*n_no_loc-8)+n_gdl
    troca = B[:,ind_B]
    B[:,ind_B] = -A[:,ind_A]
    A[:,ind_A] = -troca*GEL
    if(DESL_PR[de,5]!=0)
        n_el=DESL_PR[de,5]
        n_no_loc=DESL_PR[de,6]
        ind_B=(6*n_el+2*n_no_loc-8)+n_gdl
        A[:,ind_A]=A[:,ind_A]-GEL*B[:,ind_B]
        B[:,ind_B].=0
    end
end
	b=B*todos_valores' # C�lculo do vetor b
return A,B,DESL_PR,b
end
