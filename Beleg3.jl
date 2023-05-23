

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


for i = 1:10
    #Koeffizienten cj
    cj[i] = N[i]/A[i]
end

println(cj)

