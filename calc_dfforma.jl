function calc_dfforma(qsi)

	dN1dqsi=(2*qsi-1)/2
	dN2dqsi=-2*qsi
	dN3dqsi=(2*qsi+1)/2

	return dN1dqsi,dN2dqsi,dN3dqsi
end
