function newton(x,A,b,h,milocal,somalisttrac,nnoscontato,nlinhasA,nlinhasAA,noscontato,eps2)

contato=verifica_contato(x,h,milocal,nlinhasA,nlinhasAA,nnoscontato,somalisttrac,noscontato)
A2,b2=aplica_contato(A,b,h,milocal,somalisttrac,contato,nlinhasA,nlinhasAA,noscontato)
err=1000
x0=x
println(eps2)
while err>eps2
	println(err)
    dy=inv(A2)
    y = A2*x-b2    
    d = -dy*y # -1/derivative
    x = x+d
println(x[1:100])    
println(x0[1:100])    
    err =sqrt([x-x0]'*[x-x0])
#	aux1=[x-x0]
#	aux2=aux1'
#	err=sqrt(norm(aux2*aux1))
	   
    println("Resíduo = $err")    
    x0=x    
    contato=verifica_contato(x,h,milocal,nlinhasA,nlinhasAA,nnoscontato,somalisttrac,noscontato)
    A2,b2=aplica_contato(A,b,h,milocal,somalisttrac,contato,nlinhasA,nlinhasAA,noscontato)

end
return x,contato
end
