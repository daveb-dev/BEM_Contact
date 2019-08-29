function dad_contato()
	# Geometry input
	w=6.5
	R=70
	theta=asin(w/R)
	x=R*cos(theta)
	y=R-x
	PONTOS  =[1   -w y
		2    w y
		3    w 2.5*w
		4   -w 2.5*w
		5    w -y
		6   -w -y
		7   -w -2.5*w
		8    w -2.5*w
		]

	# Segments that defines geometry.
	# SEGMENTOS=[Node, Initial, Final, Radius, Type]
	# Radius: .> 0 -> Centre on the left of segment. 
	#         .< 0 -> Centre on the right of segment.
	#          = 0 -> Straight line.

	SEGMENTOS = [1 1 2  R
	    2 2 3  0
	    3 3 4  0
	    4 4 1  0
	    5 5 6  R
	    6 6 7  0
	    7 7 8  0
	    8 8 5  0    
	    ]

	# Mesh definition

	ref=6
	MALHA = [1  20*ref
	    2  1*ref
	    3  1*ref
	    4  1* ref
	    5  20*ref
	    6  1*ref
	    7  1*ref
	    8  1* ref]

	cargav=100; # Vertical load
	cargah=15.; # Hor load applied from left to right.
	cargah2=15.; # Hor load applied from right to left.
	# BC in segments
	# CCSeg=[Segment,BC type in x, BC value in x,BC type in y, BC value in y]
	# BC type = 0 => displacement known
	# BC type = 1 => traction known
	# Para condições de contorno de força normal conhecida proceder:
	# tipo da CDC em x = 2; valor da CDC em x = valor da força normal
	# Neste caso pode-se atribuir quaisquer valores para tipo da CDC em y e
	# para valor da CDC em y
	npassos=1

	CCSeg=[1 1 0 1 0
	    2 1 0 1 0
	    3 1 0 1 -cargav/npassos
	    4 1 0 1 0
	    5 1 0 1 0
	    6 1 0 1 0
	    7  0 0 0 0
	    8 1 0 1 0]


	nores::Int32=2*(MALHA[1,2]+MALHA[2,2]+MALHA[3,2]/2)+1
	NOS_RES=[nores 1]

	# Probable contact surfaces
	# supcontato=[seg 1, seg 2, reg 1, reg2]; 
	supcontato=[1 5]

	# Number of internal points
	npi=7
	NP_int=[npi 3*npi;npi 3*npi];  
	E=73.4e3; # Young modulus
	ni=0.33; # Poisson ratio

	mi=0.3; # Friction coefficient
	

	# Problem type
	# Type = 1 => Plane stress
	# Type = 2 => Plane strain
	tipo_problema=2


	# Output shown in colour map
	# op = "dt" => Total displacement
	# op = "sigx" => sigmax
	# op = "sigy" => sigmay
	# op = "tauxy" => tauxy
	# op = "VM" => von Mises stress

	op="VM"

	return PONTOS,SEGMENTOS,MALHA,supcontato,CCSeg,E,ni,mi,cargav,w,R,tipo_problema,NOS_RES
end
