
export afun, Pfun, Tfun, Afun, nearestto, findzero, AMfun, PMfun

#Speed of sound
"""
afun(gamma, R, T; Tconvert=false)
Calculate the speed of sound. 

### Inputs
- gamma - ratio of specific heats
- R - Specific gas constant
- T - Absolute static temperature
- Tconvert - Marks whether temperature is given in C (K if false)

### Outputs
- a - speed of sound

### Notes
Equation chapter 2? 
"""
function afun(gamma, R, T; Tconvert=false)
    if Tconvert
        T += 273.15
    end
    return sqrt(gamma*R*T)
end

#Isentropic pressure function
"""
Pfun(gamma, M)
Calculate the ratio static pressure over stagnation pressure.

### Inputs
- gamma - ratio of specific heats
- M - Mach number of the flow

### Outputs
- Pn/P0n - ratio of static pressure over stagnation pressure

### Notes
Equation 3.15
"""
function Pfun(gamma, M)
    p0_p = (1 + ((M^2)*(gamma-1)/2))^(gamma/(gamma-1))
    return 1/p0_p
end

#Isentropic temperature function
function Tfun(gamma, M)
    t0_t = (1 + ((M^2)*(gamma-1)/2))
    return 1/t0_t
end

function Afun(gamma, M)
    p1 = 2/(gamma+1)
    p2 = (1 + ((M^2)*(gamma-1)/2))
    p3 = (gamma + 1)/(2*(gamma-1))
    return ((p1*p2)^p3)/M
end

function AMfun(Gamma, Astar;subsonic=true)
    if subsonic
        rng = (0,1)
    else
        rng = (1, 200)
    end
    fun(m; gamma=Gamma, astar=Astar) = Afun(gamma, m) - astar
    mach = find_zero(fun, rng, Bisection())
    return mach
end

"""
    PMfun(Gamma, Pstar, ; subsonic=true)

### Inputs
- Gamma - ratio of specific heats of the function
- Pstar - ratio of static pressure over stagnation pressure
- subsonic - whether you'd like the subsonic solution or not

### Outputs
- M - Mach number of flow

### Notes 
Equation in chapter 2 solved numerically for M
"""
function PMfun(Gamma, Pstar;subsonic=true)
    if subsonic
        rng = (0,1)
    else
        rng = (1, 200)
    end
    fun(m; gamma=Gamma, pstar=Pstar) = Pfun(gamma, m) - pstar
    mach = find_zero(fun, rng, Bisection())
    return mach
end