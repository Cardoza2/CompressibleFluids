"""
    rfPfun(M2, γ; M1=1)

Calculates static pressure ratio between two points in Rayleigh Flow. If you'd like to use point 1 as a reference point, let M1=1. 

### Inputs
- M2::AbstractFloat - The Mach number of the second point in the flow.
- γ::AbstractFloat - The ratio of specific heats. 
- M1::AbstractFloat - The Mach number of the first point in the flow. Set to 1 automatically. 

### Outputs
- P2/P1 - The ratio of specific heats. If M1=1, then P1 is the pressure of a reference state. 

### Notes 
- Equation 10.8
"""
function rfPfun(M2, γ; M1=1)
    top = 1 +(γ*(M1^2))
    bot = 1 + (γ*(M2^2))
    return top/bot
end
export rfPfun

"""
    rfTfun(M2, γ; M1=1)

Calculates the static temperature ratio between two points in Rayleigh flow. The first point can be used as a reference point if M1 is set to one. 

### Inputs
- M2::AbstractFloat - The mach number of the second point in the flow.
- γ::AbstractFloat - The ratio of specific heats of the fluid. 
- M1::AbstractFloat - The mach number of the first (reference) point of the flow. 

### Outputs
- T2/T1 - The ratio of static temperatures. 

### Notes
- Equation 10.9
"""
function rfTfun(M2, γ; M1=1)
    top = (M2^2)*(1 +(γ*(M1^2)))^2
    bot = (M1^2)*(1 + (γ*(M2^2)))^2
    return top/bot
end
export rfTfun

function rfTMfun(TT, γ;subsonic=true)
    if subsonic
        rng = (0,1)
    else
        rng = (1, 20)
    end
    fun(m; gamma=γ, tt=TT) = rfTfun(m, gamma) - tt
    mach = find_zero(fun, rng, Bisection())
    return mach
end
export rfTMfun

"""
    rfVfun(M, γ)

Find the ratio of the velocity with respect to a reference velocity. (The book did not provide a general equation and I was too lazy to derive my own.)

### Inputs
- M::AbstractFloat - The Mach number of the fluid at the point of interest. 
- γ::AbstractFloat - The ratio of specific heats of the fluid. 

### Outputs 
- V/Vstar - The ratio of the velocity over the reference state velocity. 

### Notes
- Equation 10.13
"""
function rfVfun(M, γ)
    TTstar = rfTfun(M, γ)
    PPstar = rfPfun(M, γ)
    return TTstar/PPstar
end
export rfVfun

"""
    rfRhofun(M, γ)

Ratio of the density to a reference state density. 

### Inputs
- M::AbstractFloat - Mach number of the flow at the point of interest. 
- γ::AbstractFloat - The ratio of specific heats. 

### Outputs
- ρ/ρstar

### Notes
- Equation 10.13
"""
function rfRhofun(M, γ)
    return 1/rfVfun(M, γ)
end
export rfRhofun

"""
    rfT0fun(M, γ)

### Inputs

### Outputs

### Notes
- Equation 10.14
"""
function rfT0fun(M, γ)
    top = (1 + γ)*(M^2)*(2+((γ-1)*(M^2)))
    bot = (1 + (γ*(M^2)))^2
    return top/bot
end
export rfT0fun

function rfT0Mfun(TT, γ;subsonic=true)
    if subsonic
        rng = (0,1)
    else
        rng = (1, 50)
    end
    fun(m; gamma=γ, tt=TT) = rfT0fun(m, gamma) - tt
    mach = find_zero(fun, rng, Bisection())
    return mach
end
export rfT0Mfun

"""
    rfP0fun(M, γ)

### Inputs

### Outputs

### Notes
- Equation 10.15
"""
function rfP0fun(M, γ)
    t1 = 1 + (γ*(M^2))
    t2 = 1 + γ
    t3 = 2+((γ-1)*(M^2))
    t4 = γ/(γ-1)
    return t1*((t3/t2)^t4)/t2
end
export rfP0fun

# #This wasn't what I needed for 10.2
# function rfSfun(M, γ, TT1)
#     t1 = log(TT1)
#     t2 = (γ-1)/γ
#     t3 = 1 + (γ*(M^2))
#     t4 = 4*γ*(M^2)*TT1
#     t5 = sqrt((t3^2)-t4)
#     posi = t1 - t2*log((t3 + t5)/2)
#     negi = t1 - t2*log((t3 - t5)/2)
#     return posi, negi
# end
# export rfSfun

function rfSfun(M, M1, γ)
    t1 = (M/M1)^2
    t2t = 1 + (γ*(M1^2))
    t2b = 1 + (γ*(M^2))
    t2 = (t2t/t2b)^((γ+1)/γ)
    # println("M: ", M)
    # println("M1: ", M1)
    # println(t2)
    # println("")
    return log(t1*t2)
end
export Sfun

function Cpfun(γ, R)
    return γ*R/(γ-1)
end
export Cpfun