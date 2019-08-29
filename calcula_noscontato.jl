function calcula_noscontato(interfaces,MALHA)
	ninterfaces=size(interfaces,1)
	icont=0
	noscontato=zeros(Int64,2,interfaces[1,6])
	elemcontato=0
	for i=1:ninterfaces
		nnoscontato=interfaces[i,6] # Number of contact nodes
		nelemcontato=(nnoscontato-1)÷2
		elemcontato=zeros(Int64,2,nelemcontato)
		seg1=interfaces[i,4]
		seg2=interfaces[i,5]
		    if(seg1==1)
			eleiniseg1=1
		    else
			eleiniseg1=sum(MALHA[1:seg1-1,2])+1
		    end
	    elefimseg2=sum(MALHA[1:seg2,2])
	    noscontato[2,2*icont+1:2*icont+1+2*nelemcontato]=range(2*elefimseg2+1,2*elefimseg2-nnoscontato+2,step=-1)
	    noscontato[1,2*icont+1:2*icont+1+2*nelemcontato]=2*eleiniseg1-1:2*eleiniseg1+nnoscontato-2
	    elemcontato[2,icont+1:icont+nelemcontato]=range(elefimseg2,elefimseg2-nelemcontato+1,step=-1)
	    elemcontato[1,icont+1:icont+nelemcontato]=eleiniseg1:eleiniseg1+nelemcontato-1
	    icont=icont+nelemcontato
	end
	return noscontato, elemcontato
end
