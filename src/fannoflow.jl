

"""
    fldfun(M, γ)

Returns the value fL_max/Dh

### Inputs
- M::Float64 - Mach number of the flow
- γ::Float64 - ratio of specific heats of the flow

### Outputs
- fld - fL_max/Dh

### Notes
- Equation 9.16
- constant area, frictional, adiabatic flow
"""
function fldfun(M, γ)
    t1 = (γ+1)/(2*γ)
    t2 = log(((γ+1)/2)/(1+((γ-1)*(M^2)/2)))
    t3 = (1 - (1/(M^2)))/γ
    t4 = ((γ+1)/(2*γ))*log(1/(M^2))
    return (t1*t2) - t3 - t4
end
export fldfun

function fldMfun(fld, γ; subsonic=true)
    if subsonic
        rng = (0.000001,1) #zero is always going to be a solution... so lets avoid that solution. 
    else
        rng = (1, 200)
    end
    fun(m; gamma=γ, FLD=fld) = fldfun(m, gamma) - FLD
    mach = find_zero(fun, rng, Bisection())
    return mach
end
export fldMfun

"""
    ffPfun(M1, γ; M2=1)

Returns the static pressure ratio of the flow. 

### Inputs
- M::Float64 - Mach number of the flow
- γ::Float64 - ratio of specific heats of the flow

### Outputs
- P2/P1::Float64 - static pressure ratio

### Notes
- Equation 9.22
- Note that when setting M2 to 1, this makes this into equation 9.19 (which is equation 9.26), which references the unity state. Otherwise this function compares between two given locations.
- Note that this is backwards 
"""
function ffPfun(M1, γ; M2=1)
    t1 = M2/M1
    t2top = 2 + ((γ-1)*(M2^2))
    t2bot = 2 + ((γ-1)*(M1^2))
    t2 = sqrt(t2top/t2bot)
    println("Remember this is Pstar/P1.")
    return 1/(t1*t2)
end
export ffPfun

function ffPMfun(Pstar, Gamma;subsonic=true)
    if subsonic
        rng = (0,1)
    else
        rng = (1, 20)
    end
    fun(m; gamma=Gamma, pstar=Pstar) = ffPfun(m, gamma) - pstar
    mach = find_zero(fun, rng, Bisection())
    return mach
end
export ffPMfun

"""
    ffTfun(M1, γ;M2=1)

Returns the static temperature ratio of the flow. 

### Inputs
- M1::Float64 - Mach number of the first flow
- γ::Float64 - ratio of specific heats of the flow
- M2::Float64 - Mach number of the second flow

### Outputs
- T2/T1::Float64 - static temperature ratio

### Notes
- Equation 9.20
- Note that when setting M2 to 1, this makes this into equation 9.25, which references the unity state. Otherwise this function compares between two given locations. 
"""
function ffTfun(M1, γ;M2=1)
    top = 2 + ((γ-1)*(M2^2))
    bot = 2 + ((γ-1)*(M1^2))
    return 1/(top/bot)
end
export ffTfun

"""
    ffRhofun(M1, γ; M2=1)

Returns the static density ratio. 

### Inputs
- M1::Float64 - Mach number of the first flow
- γ::Float64 - ratio of specific heats of the flow
- M2::Float64 - Mach number of the second flow

### Outputs
- ρ2/ρ1::Float64 - static density ratio

### Notes
- Equation 9.23
- Note that when setting M2 to 1, this makes this into equation 9.27, which references the unity state. Otherwise this function compares between two given locations. 
"""
function ffRhofun(M1, γ; M2=1)
    t2 = ffPfun(M1, γ; M2=M2)
    t1 = ffTfun(M1, γ; M2=M2)
    return 1/(t1/t2)
end
export ffRhofun

"""
    ffP0fun(M1, γ; M2=1)

Returns the stagnation pressure ratio. 

### Inputs
- M1::Float64 - Mach number of the first flow
- γ::Float64 - ratio of specific heats of the flow
- M2::Float64 - Mach number of the second flow

### Outputs
- P02/P01::Float64 - stagnation pressure ratio

### Notes
- Equation 9.24
- Note that when setting M2 to 1, this makes this into equation 9.28, which references the unity state. Otherwise this function compares between two given locations. 
"""
function ffP0fun(M1, γ; M2=1)
    t1 = (M2/M1)
    t2top = 2 + ((γ-1)*(M1^2))
    t2bot = 2 + ((γ-1)*(M2^2))
    t2 = (t2top/t2top)^((γ+1)/(2*(γ-1)))
    return 1/(t1*t2)
end
export ffP0fun

"""
    function ffSfun(M1, γ; M2=1, R=287)

Returns the change in entropy between two frictional flows. 

### Inputs
- M1::Float64 - Mach number of the first flow
- γ::Float64 - ratio of specific heats of the flow
- M2::Float64 - Mach number of the second flow

### Outputs
- Δs::Float64 - the change in entropy in the flow between two points (s2-s1). 

### Notes
- Equation 4.14
- Note that when setting M2 to 1 makes this function compute the change in entropy between the reference state and the current point. 
"""
function ffSfun(M1, γ; M2=1, R=287)
    return -R*(log(ffP0fun(M1, γ; M2=M2)))
end

"""
    ffVfun(M1, γ; M2=1)

Returns the ratio of the fluid velocities. 

### Inputs
- M1::Float64 - Mach number of the first flow
- γ::Float64 - ratio of specific heats of the flow
- M2::Float64 - Mach number of the second flow

### Outputs
- V2/V1::Float64 - the change in entropy in the flow between two points (s2-s1). 

### Notes
- Equation 9.23
- Note that when setting M2 to 1 makes this function compute the change in entropy between the reference state and the current point. 
"""
function ffVfun(M1, γ; M2=1)
    return 1/ffRhofun(M1, γ; M2=M2)
end
export ffVfun  

function fun933(M1; γ=1.4)
    t1 = (γ+1)/γ
    t2top = 2 + ((γ-1)*(M1^2))
    t2bot = (γ+1)*(M1^2)
    t2 = log(t2top/t2bot)
    t3top = 2*(1+(γ*(M1^2)))*((M1^2)-1)
    t3bot = γ*(M1^2)*(2+((γ-1)*(M1^2)))
    t3 = t3top/t3bot
    return t1*t2 + t3
end
export fun933
"""
    fun952(M, A, dAdx, D, f, γ)

Calculates the derivative of the Mach number with respect to x. 

### Inputs
- M - Mach number of location
- A - Area of location
- dAdx - derivative of how the area is changing 
- D - Hydraulic Diameter
- f - Fanning friction factor
- γ - ratio of specific heats

### Notes
- Equation 9.52
- Assuming dAdx is constant
"""
function fun952(M, A, dAdx, D, f, γ)
    t1 = (((γ-1)*(M^2))+2)*M/2
    t2top = -(dAdx/A) + γ*(M^2)*f/(2*D)
    t2bot = 1 - (M^2)
    t2 = t2top/t2bot
    return t1*t2
end
export fun952

function ChangingAreaDuct()
end
export ChangingAreaDuct


