function format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg,E,ni,tipo_problema)

# Programa para formatacao dos dados de entrada
#
#   Autor: Eder Lima de Albuquerque


cont_nos = 0  # Contador do n�mero de n�s
cont_el = 0	# Contador do n�mero de elementos
num_seg= length(SEGMENTOS[:,1])	# N�mero de segmentos nos contornos
p_ini = SEGMENTOS[1,2]
NOS=zeros(552,3)
ELEM=zeros(276,4)
CDC=zeros(276,13)
#NOS=[]
#ELEM=[]
#CDC=[]

dNdqsi=zeros(3,3)

#______________________________________________________________________
# Defini��o da maior dimens�o do problema
max_dl = 0
for lin = 1 : num_seg
    p1 = SEGMENTOS[lin,2]
    p2 = SEGMENTOS[lin,3]
    xp1 = PONTOS[p1,2];	yp1 = PONTOS[p1,3]
    xp2 = PONTOS[p2,2];	yp2 = PONTOS[p2,3]
    dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2)
    if dl .> max_dl
        max_dl = dl
    end
end
#_____________________________________________________________________

no_ini=1
t=1
p2=0
no1_prox=0
while(t<num_seg)  	# While sobre todos os segmentos
    while(p2!=p_ini)
        num_el_lin = MALHA[t,2]	# N�mero de elementos no segmento t
        
        
        # Coordenadas iniciais e finais dos pontos de cada segmento
        # (x1l,y1l,x2l,y2l)
        p1  = SEGMENTOS[t,2]
        p2  = SEGMENTOS[t,3]
        x1l = PONTOS[p1,2]
        y1l = PONTOS[p1,3]
        x2l = PONTOS[p2,2]
        y2l = PONTOS[p2,3]
                
        # 1. Gera��o das matrizes NOS; NOS_GEO; NOS_DRM; ELEM e ELEM_GEO
        if(SEGMENTOS[t,4]==0) # O segmento � uma linha reta
            # Incremento nas dire��es x e y
            delta_x = x2l - x1l
            delta_y = y2l - y1l
        else # O segmento � um arco
            # C�lcula as coordenadas do centro do arco
            r = SEGMENTOS[t,4]
            xc,yc=calcula_centro(x1l,y1l,x2l,y2l,r)
            # Dist�ncia entre p1 e c[r1] e entre p2 e c[r2]
            r1 = sqrt((x1l-xc)^2+(y1l-yc)^2)
            r2 = sqrt((x2l-xc)^2+(y2l-yc)^2)
            if abs(r1-r2)<.00001*max_dl
                # C�lcula o �ngulo entre as linhas do ponto c at� p1[tet1]
                #  e entre c at� p2[tet2]
                tet1,tet2 = calcula_arco(x1l,y1l,x2l,y2l,xc,yc,r)
                if tet2 < tet1
                    tet2 = tet2 + 2*pi
                end
                
                # �ngulo entre os setores definidos pelo arco
                if SEGMENTOS[t,4] .> 0
                    tet = abs(tet2-tet1)
                    sig = 1
                else
                    tet = 2*pi-abs(tet2-tet1)
                    sig = -1
                end
                
                # �ngulo entre dois n�s
                divtet = tet/(2*num_el_lin)
            else
                error("Error in the data input file: Wrong central point")
            end
        end
        
        
        # Gera��o dos elementos e n�s
        for i = 1 : num_el_lin
            
            if(SEGMENTOS[t,4]==0) # O segmento � uma linha reta
                x_i = x1l + delta_x/num_el_lin*(i-1)			# Coordenada inicial x do elemento
                y_i = y1l + delta_y/num_el_lin*(i-1)			# Coordenada inicial y do elemento
                x_m = x1l + delta_x/num_el_lin*(i-.5)	# Coordenada x do ponto central do elemento
                y_m = y1l + delta_y/num_el_lin*(i-.5)	# Coordenada y do ponto central do elemento
                lx=x_m-x_i                              # Dist�ncia na dire��o x entre os n�s 1 e 2
                ly=y_m-y_i                              # Dist�ncia na dire��o y entre os n�s 1 e 2
                
            else  # O segmento � um arco
                # Calcula as coordenadas dos n�s
                x_i = xc+r1*cos(tet1+2*(i-1)*sig*divtet)
                y_i = yc+r1*sin(tet1+2*(i-1)*sig*divtet)
                x_m = xc+r1*cos(tet1+(2*i-1)*sig*divtet)
                y_m = yc+r1*sin(tet1+(2*i-1)*sig*divtet)
            end
            cont_el = cont_el + 1
            # Atribui coordenadas aos n�s
            if(no1_prox==0) # O primeiro do elemento n� precisa ser criado
                cont_nos = cont_nos + 1
                NOS[cont_nos,:]=[cont_nos,x_i,y_i]
                no1=cont_nos
            else
                no1=no1_prox
            end
            cont_nos = cont_nos + 1
            NOS[cont_nos,:]=[cont_nos,x_m,y_m]
            no2=cont_nos
            if(p2!=p_ini || i<num_el_lin)
                cont_nos = cont_nos + 1
                if(SEGMENTOS[t,4]==0) # Linha reta
                    x_f=x_m+lx
                    y_f=y_m+ly
                else                 # Arco
                    x_f = xc+r1*cos(tet1+(2*i)*sig*divtet)
                    y_f = yc+r1*sin(tet1+(2*i)*sig*divtet)
                end
                NOS[cont_nos,:]=[cont_nos,x_f,y_f]
                no3=cont_nos
                no1_prox=no3
            else
                if(no_ini==0)
                    cont_nos = cont_nos + 1
                    if(SEGMENTOS[t,4]==0) # Linha reta
                        x_f=x_m+lx
                        y_f=y_m+ly
                    else                 # Arco
                        x_f = xc+r1*cos(tet1+(2*i)*sig*divtet)
                        y_f = yc+r1*sin(tet1+(2*i)*sig*divtet)
                    end
                    NOS[cont_nos,:]=[cont_nos,x_f,y_f]
                    no3=cont_nos
                    no1_prox=0
                else
                    no3=no_ini
                    no1_prox=0
                end
            end
            ELEM[cont_el,:]=[cont_el,no1,no2,no3]
        end # fim do for i = 1 : num_el_lin

        if p2 == p_ini
            if t < num_seg
                p_ini = SEGMENTOS[t+1,2]
                if(SEGMENTOS[t+1,3]==2)
                    no_ini = 0
                else
                    no_ini = cont_nos+1
                end
            end
        end
        t=t+1
    end                                  # fim do while p2
end                                   # fim do while(t<num_lin)
#

# Gera��o da matriz CDC [Condi��es de Contorno]
# CDC = [n. do elemento, tipo de cdc, valor da cdc]
# Tipos de cdc: 0 : O deslocamento � conhecido
#               1 : A for�a de superf�cie � conhecida
cont_el2 = 0
qsi_nos=[-1,0,1]
for i=1:3
    dN1dqsi,dN2dqsi,dN3dqsi=calc_dfforma(qsi_nos[i]) # Calcula as derivadas das
    dNdqsi[i,1:3]=[dN1dqsi,dN2dqsi,dN3dqsi]
    #    fun��es de forma
end

for l = 1 : length(SEGMENTOS[:,1])
    n_el_seg = MALHA[l,2]
    el_ini = cont_el2 + 1
    el_fin = cont_el2 + n_el_seg
    tipoCDCx=CCSeg[l,2]
    if(tipoCDCx==2)
        ForcaNormal=1
        valor=CCSeg[l,3]
    else
        ForcaNormal=0
    end
    valorCDCx=[CCSeg[l,3] CCSeg[l,3] CCSeg[l,3]]
    tipoCDCy=CCSeg[l,4]
    valorCDCy=[CCSeg[l,5] CCSeg[l,5] CCSeg[l,5]]
    for el = el_ini : el_fin
        if(ForcaNormal==1)
            no1 = ELEM[el,2]
            no2 = ELEM[el,3]
            no3 = ELEM[el,4]
            x1 = NOS[no1,2];	y1 = NOS[no1,3]
            x2 = NOS[no2,2];	y2 = NOS[no2,3]
            x3 = NOS[no3,2];	y3 = NOS[no3,3]
            tipoCDCx=1
            tipoCDCy=1
            for i=1:3
                dxdqsi= dNdqsi[i,1]*x1+dNdqsi[i,2]*x2+dNdqsi[i,3]*x3
                dydqsi= dNdqsi[i,1]*y1+dNdqsi[i,2]*y2+dNdqsi[i,3]*y3
                dgamadqsi=sqrt(dxdqsi^2+dydqsi^2)
                sx=dxdqsi/dgamadqsi; # Componente x do vetor tangente na extremidade do elemento
                sy=dydqsi/dgamadqsi; # Componente y do vetor tangente na extremidade do elemento
                n1=sy; # Componente x do vetor normal unit�rio na extremidade do elemento
                n2=-sx; # Componente y do vetor normal unit�rio na extremidade do elemento
                valorCDCx[i]=valor*n1
                valorCDCy[i]=valor*n2
            end
        end
        CDC[el,:] = [el,tipoCDCx,valorCDCx[1],tipoCDCy,valorCDCy[1],      tipoCDCx,valorCDCx[2],tipoCDCy,valorCDCy[2],tipoCDCx,valorCDCx[3],tipoCDCy,valorCDCy[3]]
	    end
	    cont_el2 = el_fin
	end

	if(tipo_problema==1) # Se tipo = 1 => Estado plano de tens�o
	    E=E/(1-ni^2)
	    ni=ni/(1-ni)
	end
	return NOS,ELEM,CDC,E,ni
end
