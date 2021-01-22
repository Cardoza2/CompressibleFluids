
"""
## PMangle(M, γ)
Calculates the Prandtl-Meyer angle of the flow

### Inputs
- M - Flow Mach number
- γ - fluid ratio of specific heats

### Outputs
- ν - Prandtl-Meyer angle

### Notes
Equation 7.10, assumes isentropic flow.
"""
function PMangle(M, γ)
    gp = γ+1
    gn = γ-1
    t1 = sqrt(gp/gn)
    t2 = sqrt(gn*((M^2)-1)/gp)
    t3 = sqrt((M^2)-1)
    return t1*atand(t2) - atand(t3)
end
export PMangle

"""
#### efPfun(M1, M2, γ)
Calculates the static pressure ratio before and after an expansion fan.

### Inputs
- M1 - Mach number of incoming flow
- M2 - Mach number of exiting flow
- γ - ratio of specific heats of the fluid

### Outputs
- P2/P1 - ratio of static pressures

### Notes
Equation 7.13
"""
function efPfun(M1, M2, γ)
    top = (γ-1)*(M1^2)/2 + 1
    bot = (γ-1)*(M2^2)/2 + 1
    return (top/bot)^(γ/(γ-1))
end
export efPfun

#TODO: add a efPMfun that is just PMfun from isentropic flow

"""
#### efTfun(M1, M2, γ)
Calculates the ratio of static Temperatures based on the Mach numbers.

### Inputs
- M1 - Mach number of incoming flow
- M2 - Mach number of exit flow
- γ - ratio of specific heats for fluid

### Outputs
- T2/T1 - ratio of static temperatures 

### Notes
Equation 7.14, Note that equations 7.13 (efPfun) and 7.14 vary by an exponent.
"""
function efTfun(M1, M2, γ)
    #Validated
    top = (γ-1)*(M1^2)/2 + 1
    bot = (γ-1)*(M2^2)/2 + 1
    return (top/bot)
end
export efTfun




"""
#### MPMfun(ν, γ)
Calculate the Mach number from the PM angle

### Inputs
- ν - Prandtl-Meyer Angle
- γ - ratio of specific heats

### Outputs
- M - Mach number 

### Notes
Equation 7.10 solved for M using a root finder, assuming isentropic flow
"""
function MPMfun(ν, γ)
    rng = (1,100)
    fun(m; gamma=γ, nu=ν) = PMangle(m, gamma) - nu
    mach = find_zero(fun, rng, Bisection())
    return mach
end
export MPMfun

"""
#### expansionfan(M1, Δ, γ ; P1=0, T01=0)
Compute the aspects of an expansion fan

### Inputs
- M1 - Incoming flow mach number
- Δ - Turning angle. Angle "below the horizon".
- γ -  ratio of specific heats

### Outputs
A fan dictionary containing all of the input and output values, including:
- Mn - Mach number
- Δ - Turning angle
- γ - ratio of specific heats
- νn - Prandtl-Meyer angle
- P2P01 - ratio of static pressure over stagnation pressure
- T2T01 - ratio of static temperature over stagnation temperature
- P2 - static pressure after the expansion fan
- T2 - static temperature after the expansion fan

### Notes
using equation 7.10 and PM angles mainly to solve, and some isentropic flow relations from chapter 2 to provide further information. 
"""
function expansionfan(M1, Δ, γ; P1=0, T1=0)
    fan = Dict()
    fan["M1"] = M1
    fan["Δ"] = Δ
    fan["γ"] = γ
    fan["ν1"] = PMangle(M1, γ)
    fan["ν2"] = Δ+fan["ν1"]
    fan["M2"] = MPMfun(fan["ν2"], γ)
    fan["P2P01"] = Pfun(γ, fan["M2"])
    fan["T2T01"] = Tfun(γ, fan["M2"])
    fan["P2P1"] = efPfun(M1, fan["M2"], γ)
    fan["T2T1"] = efTfun(M1, fan["M2"], γ)
    #Note that I might need N2N01 or N2N1 in order to solve for N2(N for P or T).

    if P1>0
        fan["P2"] = fan["P2P1"]*P1
    end
    if T1>0
        fan["T2"] = fan["T2T1"]*T1
        fan["T1"] = T1
    end
    return fan
end
export expansionfan



"""
#### efPMfun(Prat, M1, γ)
Calculates the Mach number after an expansion fan based on the static pressure ratio before and after an expansion fan, and the mach number of the incident flow.

### Inputs
- Prat - Static pressure ratio
- M1 - Incoming flow Mach number
- γ - ratio of specific heats

### Outputs
- M - Mach number of flow

### Notes 
Equation 7.13 solved numerically for M2
"""
function efPMfun(Prat, M1, γ)
    rng = (1, 100)
    fun(m2; m1=M1, gamma=γ, prat=Prat) = efPfun(m1, m2, gamma) - prat
    mach2 = find_zero(fun, rng, Bisection())
    return mach2
end
export efPMfun

