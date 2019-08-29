function foo()

global x0
global somalisttrac
milocal=mi
NOSt[:,2]=NOS[:,2]+desl[:,2] 
NOSt[:,3]=NOS[:,3]+desl[:,3] 
h=calc_dist(NOSt,noscontato)
println("Newton")
erromax = 2e-10
x,contato=newton(x0,A2,b2,h,milocal,somalisttrac,nnoscontato,nlinhasA,nlinhasAA,noscontato,erromax)

end # End of function foo
