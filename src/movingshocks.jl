export msRhofun, msPfun, msTfun, msVSfun, msSfun

"""
msRhofun(s, v; w=0)

Inputs:
- s - Speed of the shock
- v - speed of the flow after the shock
- w - speed of the flow which the shock is propagating into

Outputs: 
- rho1/rho2 - Density ratio *** NOTE THIS IS BACKWARDS FROM OTHERS ***

notes:
    Equation 5.6 . For static properties -> works for a moving or stationary frame of reference. 
"""
function msRhofun(s, v; w=0)
    S = s-w
    V = v-w
    return 1 - V/S
end

"""
mspfun(s, v, a1, γ; w=0)

inputs
    s - Speed of the shock
    v - speed of the flow after the shock
    a1 - speed of sound in which the shock is propagating into
    γ - ratio of specific heats of the gas
    w - speed of the fluid of the gas which the shock propagates into
returns 
    p2/p1 - ratio of static pressures before and after the wave

notes
    Equation 5.7
"""
function msPfun(s, v, a1, γ; w=0)
    S = s-w
    V = v-w
    t2 = γ*S*V/(a1^2)
    return 1 + t2
end

"""
mspfun(s, a1, γ; w=0)

inputs
    s - Speed of the shock
    a1 - speed of sound in which the shock is propagating into
    γ - ratio of specific heats of the gas
    w - speed of the fluid of the gas which the shock propagates into
returns 
    p2/p1 - ratio of static pressures before and after the wave

notes
    Equation 5.8
"""
function msPfun(s, a1, γ; w=0)
    S = s-w
    t1 = 2*γ/(γ+1)
    t2 = (S/a1)^2
    t3 = (γ-1)/(γ+1)
    return (t1*t2)-t3
end

"""
mspfun(v, a1; γ)

inputs
    a1 - speed of sound in which the shock is propagating into
    γ - ratio of specific heats of the gas
    w - speed of the fluid of the gas which the shock propagates into
returns 
    p2/p1 - ratio of static pressures before and after the wave

notes
    Equation 5.12
"""
function msPfun(v, a1; γ=1.4, w=0)
    V = v-w
    t2 = γ*((γ+1)/4)
    t3 = (V/a1)^2
    t4 = γ*V/a1
    t5_1 = ((γ+1)/4)^2
    t5_2 = (V/a1)^2
    t5 = sqrt(1 + (t5_1*t5_2))
    return 1 + t2*t3 + t4*t5
end

"""
mstfun(a1, a2)

inputs
    a1 - speed of sound of incident gas
    a2 - speed of sound of following gas
returns
    T2/T1 - ratio of temperatures

notes
    equation 5.13
"""
function msTfun(a1,a2)
    return (a2^2)/(a1^2)
end

"""
mstfun(s, v, a1, γ; w=0)

inputs
    s - speed of shock
    v - speed of proceeding fluid
    a1 - speed of sound of incident gas
    γ - ratio of specific heats
    w = speed of incident gas

returns 
    T2/T1 - ratio of temperatures across shock

notes
    equation 5.13
"""
function msTfun(s, v, a1, γ; w=0)
    S = s-w
    V = v-w
    t2 = (γ-1)/(a1^2)
    t3 = S*V - (V^2)/2
    return 1 + t2*t3
end

"""
msvsfun(s, a1, γ; w=0)

Inputs:
- s - speed of shock
- a1 - speed of sound of incident gas
- γ - ratio of specific heats 
- w - speed of incident gas (positive in same direction as shock)
Outputs: 
- V/S - ratio of velocity of proceeding gas to shock speed

notes
    equation 5.9, This takes in and returns the s, and v that aren't adjusted to have a moving reference frame. 
"""
function msVSfun(s, a1, γ; w=0)
    S = s-w
    t1 = 2/(γ+1)
    t3 = (a1/S)^2
    V_S = t1*(1-t3)
    V = V_S*S #Solve for V
    v = V + w #Transition back to V without w taken into account
    return v/s
end

"""
msSfun(v, a1, γ; w=0)

Inputs:
- v - velocity of the fluid proceeding the shock
- a1 - speed of sound of incident gas
- γ - ratio of specific heats
- w - velocity of incident gas (positive in direction of shock travel)

Outputs: 
- S - speed of shock

Notes:
    Equation 5.10
"""
function msSfun(v, a1, γ;w=0)
    V = v-w
    t1 = ((γ+1)/4)*V
    t2 = t1^2
    t3 = a1^2
    S =  t1 + sqrt(t2+t3)
    s = S + w #Shift back to moving reference frame
    return s
end

"""
msSfun(a1, P2P1; γ1 = 1.4)

Inputs:
- a1 - speed of sound of incident gass
- γ - ratio of specific heats of incident gas
- P2P1 - ratio of pressures across shockwave
- w - speed of the incident gas (positive in direction of shock)

Output:
- S - Speed of sound of the shockwave

Notes:
    Equation 5.32 (Equation 5.8 solved for S)
"""
function msSfun(a1, P2P1; γ=1.4, w=0)
    t2 = (γ+1)/(2*γ)
    t4 = (γ-1)/(2*γ)
    S = a1*sqrt(t2*P2P1 + t4)
    return S+w
end

"""
msT02T01fun(T2T1, v, a2, γ)

Inputs:
- T2T1 - ratio of static temperatures, Equation 5.13
- v - velocity of following gas
- a2 - speed of sound of following gas
- γ - ratio of specific heats
Returns:
- T02T01 - Ratio of stagnation temperatures

Notes:
Derived using equation 5.13 and ratios
"""
function msT02T01fun(T2T1, v, a2, γ;w=0)
    V = v-2
    t3 = (γ-1)/2
    t4 = (V/a2)^2
    return T2T1*(1 + t3*t4)
end
export msT02T01fun

"""
msT02T01(s, T1, γ, R; w=0)

Inputs:
- s - speed of the shock
- T1 - Temperature of incident gas
- γ - ratio of specific heats
- R - gas specific gas constant
- w - velocity of incident gas
Outputs: 
- T02/T01 - ratio of stagnation temperatures
"""
function msT02T01(s, T1, γ, R; w=0)
    S = s-w
    a1 = afun(γ, R, T1)
    V_S = msVSfun(S, a1, γ)
    V = V_S*S #I'm not 100% on this line, guess I ought to test it. 
    T2T1 = msTfun(S, V, a1, γ)
    T2 = T2T1*T1
    a2 = afun(γ, R, T2)
    return msT02T01fun(T2T1, V, a2, γ)
end
export msT02T01