using Distributions
using QuadGK

#Erstellen Arrays für die Stäbe

#Wanddicke
t = [6.3 6.3 6.3 6.3 6.3 8.0 6.3 6.3 6.3 10.0]
#Querschnittsfläche
A = [3.77 3.77 3.77 3.77 3.77 4.70 3.77 3.77 3.77 5.74] * 10^-3
#Flächenträgheitsmoment
I = [1.46 1.46 1.46 1.46 1.46 1.78 1.46 1.46 1.46 2.10] * 10^-5
#Länge
l = [5.6 4.2 4.2 7.0 4.2 sqrt(35.28) 5.6 4.2 7.0 4.2]

#Stabkräfte (V=1)
N = [0 1 0 1.25 -0.75 -1.414 -1 0.75 0 -1.175]
E = 2.1e8

#Belastung
muh_x1 = 410
sigma_x1 = 70

#Material
muh_x2 = 30.2e4
sigma_x2 = 2.44e4
x02 = 19.9e4

#leere Arrays
cj = zeros(10)
x01 = zeros(10)
Pf = zeros(10)
f1_array = zeros(10)
F_min_array = zeros(10)
untere_Grenze = zeros(10)
Pf_ges  = 0
#Koeffizienten cj
for i = 1:10
    cj[i] = abs(N[i]/A[i])
end
cj = - cj
println("cj: ",cj)
println()


#Werte a und b
a = 1/sigma_x1*pi/sqrt(6)
b = muh_x1 - 0.5571/a

println("a: ",a,"  b: ", b)



#Funktionen f1_x1 und F_min_x1

#belegen des arrays x01
for i = 1:10
    if cj[i] == 0
        x01[i] = 0
    else
        x01[i] = -x02/cj[i]
    end
end
println("x01: ", x01)
println()

#Berechnung von sigma_u und muh_u
sigma_u = sqrt(log(1+(sigma_x2/(muh_x2-x02))^2))
muh_u = log(muh_x2-x02)-sigma_u^2/2
println("sigma_u :",sigma_u,"   muh_u: ",muh_u)
println()


function f_x1(x)
    return f_x1 = a*exp(-a*(x-b)*exp(-a*(x-b)))
end

function F_min_x1(x)
        z = (log(x-x02)-muh_x2)/sigma_x2
        return pdf(Normal(), z)
end

for i = 1:10

   if cj[i] == -0.0

   else
       Pf[i], error = quadgk(x -> f_x1(x) * F_min_x1(-cj[i]*x), -x02/cj[i] , 1234)
       global Pf_ges = Pf_ges + Pf[i]
   end


end

println(Pf_ges)
# result, error = quadgk(x -> f_x1(x) * F_min_x1(x), , Inf)
# println("Ergebnis der Integration: ", result)
# println("Fehler der Integration: ", error)
