function calc_solfund(x,y,x0,y0,nx,ny,E,ni) 

#Calcula as solu��es fundamentais

	GEL=E÷(2*(1+ni))
	R=sqrt((x-x0)^2+(y-y0)^2) # Raio (dist�ncia entre ponto fonte e 
		                   # ponto campo)
	Rx=(x-x0) # Componente x do raio
	Ry=(y-y0) # Componente y do raio
	# C�lculo de fatores comuns que se repetem nas solu��es

	rd1=Rx/R
	rd2=Ry/R
	prod1 = 4*pi*(1-ni)
	prod2 = (3-4*ni)*log(1/R)
	prod3 = rd1*nx + rd2*ny

	# Solu��o fundamental de deslocamento
	u11 = (prod2 + rd1^2)/(2*prod1*GEL)
	u22 = (prod2 + rd2^2)/(2*prod1*GEL)
	u12 = (rd1*rd2)/(2*prod1*GEL)

	# Solu��o fundamental de tra��o
	p11 = -(prod3*((1-2*ni)+2*rd1^2))/(prod1*R)
	p22 = -(prod3*((1-2*ni)+2*rd2^2))/(prod1*R)
	p12 = -((prod3*2*rd1*rd2)-(1-2*ni)*(rd1*ny-rd2*nx))/(prod1*R)
	p21 = -((prod3*2*rd1*rd2)-(1-2*ni)*(rd2*nx-rd1*ny))/(prod1*R)

	# Atribui��o de uij e pij �s vari�veis de sa�da u_est e p_est
	uast = [u11 u12;u12 u22]

	tast = [p11 p12;p21 p22]

	return uast,tast
end
