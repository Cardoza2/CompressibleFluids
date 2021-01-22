export Helium20, He, Turpentine, Lead, Air, Water, O2

"""
Dictionaries of different gasses and their properties. Compressiblity functions.
Adam Cardoza
Fall 2020
"""

Helium20 = Dict()
Helium20["γ"] = 1.667 #Unitless, 20°C
Helium20["R"] = 2076.9 #J/kg⋅K
Helium20["a"] = 1015 #m/s
Helium20["TK"] = 293.15 #K
Helium20["k_s"] = 5919e-9 #1/Pa
Helium20["Β_s"] = 169.0e3 #Pa
Helium20["ρ"] = 0.16 #kg/m^3

He = Dict()
He["γ"] = 5/3 #Unitless, 20°C
He["R"] = 2076.9 #J/kg⋅K
He["Β_s"] = 169.0e3 #Pa
He["Mass"] = 28.97

Turpentine = Dict() #Liquid!!!!
Turpentine["ρ"] = 870 #kg/m^3
Turpentine["k_s"] = 0.736e-9 #1/Pa

Lead = Dict() #Solid!!!
Lead["ρ"] = 11300 #kg/m^3
Lead["k_s"] = 0.061e-9 #1/Pa
Lead["Β_s"] = 16.27e9 #Pa

Air = Dict()
Air["γ"] = 1.402
Air["R"] = 287

Water = Dict()
Water["γ"] = 1.004
Water["ρ"] = 998
Water["k_s"] = 0.457e-9
Water["Β_s"] = 2.19e9 

O2 = Dict()
O2["γ"] = 1.4
O2["R"] = 0.2598e3

H = Dict()
H["γ"] = 1.4
H["R"] = 4124

