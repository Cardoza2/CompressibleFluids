"""
osPfun(M1, θ, γ; M2=false, δ=0)

Inputs: 
- M1 - Mach number of incoming flow
- θ - Oblique shock wave angle
- γ - ratio of specific heats
- M2 - If the calculation is for M2
- δ - turning angle 

Outputs:
- P2P1 - ratio of static pressures

Notes:
- Equation 6.10
"""
function osPfun(M1, θ, γ; M2=false, δ=0)
    if M2
        theta = θ-δ
        M = M2
    else
        theta = θ
        M = M1
    end

    t1 = 2*γ*(M^2)*(sind(theta)^2)
    t2 = γ+1
    t3 = γ-1
    return (t1-t3)/t2
end
export osPfun

"""
osRhofun()
equation 6.11
"""

"""
osTfun(M1, θ, γ; M2=false, δ=0)

Inputs: 
- M1 - Mach number of incoming flow
- θ - Oblique shock wave angle
- γ - ratio of specific heats
- M2 - If the calculation is for M2
- δ - turning angle 

Outputs:
- T2T1 - ratio of static temperatures

Notes:
Equation 6.12
"""
function osTfun(M1, θ, γ; M2=false, δ=0)
    if M2
        theta = θ-δ
        M = M2
    else
        theta = θ
        M = M1
    end
    t1 = (γ-1)*(M^2)*(sind(theta)^2)/2 + 1
    t2 = 2*γ*(M^2)*(sind(theta)^2)/(γ-1) - 1
    t3 = ((γ+1)^2)*(M^2)*(sind(theta)^2)/(2*(γ-1))
    return t1*t2/t3
end
export osTfun


"""
deltafun(θ, M1, γ)

Inputs:
- θ - Oblique Shock angle (defined from horizontal)
- M1 - Mach number of incoming flow
- γ - ratio of specific heats for fluid

Outputs:
- tanδ - the tangent of the turning angle

Notes: 
Equation 6.18
"""
function deltafun(θ, M1, γ)
    top = ((M1^2)*(sind(θ)^2))-1
    t1 = (M1^2)*(γ+1)/2
    t2 = ((M1^2)*(sind(θ)^2))-1
    bot = t1-t2
    tanδ = cotd(θ)*top/bot
    return tanδ
end
export deltafun

"""
rearrangeddeltafun(θ, M1, γ, δ)
equation 6.18 but rearranged. 
"""
function rearrangeddeltafun(θ, M1, γ, δ)
    theta = θ*pi/180
    delta = δ*pi/180
    x = cot(theta)
    A = (M1^2)-1
    B = ((γ+1)/2)*(M1^4)*tan(delta)
    term = (M1^2)*(γ+1)/2
    C = (1 + term)*tan(delta)
    return x^3 + C*(x^2) - A*x + (B-(A*C))
end



"""
collarsmethod(δ, M1, γ)

Inputs: 
- δ - the turning angle
- M1 - Mach number of incoming flow
- γ - ratio of specific heats of fluid. 

Outputs:
- thetas - The two solutions allowed for a θ - δ relationship. The first solution is the weak shock solution. The second solution is the strong shock solution. 
"""
function collarsmethod(δ, M1, γ)
    thetas = collect(δ:.001:90)
    y = zeros(length(thetas))
    for i = 1:length(thetas)
        y[i] = rearrangeddeltafun(thetas[i], M1, γ, δ)
    end
    ya = vcat(y, [0])
    yb = vcat([0], y)
    yc = zeros(length(ya))
    idxs = []
    for i = 1:length(ya)
        yc[i] = ya[i]*yb[i]
        if yc[i]<0
            push!(idxs, i-1)
        end
    end
    # println(length(idxs))
    roots = zeros(2)
    fun(theta; m1=M1, gamma=γ, delta=δ) = rearrangeddeltafun(theta, m1, gamma, delta)
    for i = 1:2
        # println((thetas[idxs[i]],thetas[idxs[i]+1]))
        roots[i] = find_zero(fun, (thetas[idxs[i]],thetas[idxs[i]+1]), Bisection())
    end
    return roots
end
export collarsmethod

"""
thetamax(M1, γ)

Inputs:
- M1 - Inflow mach number
- γ - ratio of specific heats

Outputs:
- θmax - The maximum theta before the flow separates 

Notes:
- Equation 6.24
"""
function thetamax(M1, γ)
    t1 = 1/(γ*(M1^2))
    t2 = (M1^2)*(γ+1)/4
    t3 = (M1^4)*(γ+1)/16
    t4 = (M1^2)*(γ-1)/2
    t5 = sqrt((γ+1)*(t3+t4+1))
    eq = t1*(t2 - 1 + t5)
    return asind(sqrt(eq))
end
export thetamax

"""
obliqueshock(delta, M1, γ)

"""
function weakobliqueshock(delta, M1, gamma; T1=0, R=0, P1=0)
    shock = Dict()
    shock["M1"] = M1
    shock["thetas"] = collarsmethod(delta, M1, gamma)
    theta = shock["thetas"][1]
    shock["M1n"] = M1*sind(theta)
    shockn = normshock(gamma, shock["M1n"])
    shock["M2n"] = shockn["M2"]
    shock["P2P1"] = shockn["P2P1"]
    shock["T2T1"] = shockn["T2T1"]
    shock["P02P01"] = shockn["P02P01"]
    shock["M2"] = shock["M2n"]/(sind(theta-delta))
    shock["M1t"] = M1*cosd(theta)
    shock["M2t"] = shock["M2"]*cosd(theta-delta)
    shock["θmax"] = thetamax(M1, gamma)
    tanδmax = deltafun(shock["θmax"], M1, gamma)
    shock["δmax"] = atand(tanδmax)

    # shock["T2T1"] = osTfun(M1, theta, gamma)
    # shock["P2P1"] = osPfun(M1, theta, gamma)
    shock["staticT"] = false
    if (T1>0 && R>0)
        # println("ran this")
        shock["staticT"] = true
        shock["T1"] = T1
        shock["T2"] = shock["T2T1"]*T1
        shock["a1"] = afun(gamma, R, shock["T1"])
        shock["a2"] = afun(gamma, R, shock["T2"])
        shock["V1"] = M1*shock["a1"]
        shock["V2"] = shock["M2"]*shock["a2"]
        shock["V1n"] = shock["M1n"]*shock["a1"]
        shock["V1t"] = shock["M1t"]*shock["a1"] 
        shock["V2n"] = shock["M2n"]*shock["a2"]
        shock["V2t"] = shock["M2t"]*shock["a2"]
    end
    shock["staticP"] = false
    if P1>0
        shock["P1"] = P1
        shock["staticP"] = true
        shock["P2"] = shock["P2P1"]*P1
    end

    return shock
end
export weakobliqueshock

function printshock(shock)
    println("Side 1:")
    println("   M: ", shock["M1"])
    println("   M1n: ", shock["M1n"])
    println("   M1t: ", shock["M1t"])
    if shock["staticT"]
        println("   V1: ", shock["V1"])
        println("   V1n: ", shock["V1n"])
        println("   V1t: ", shock["V1t"])
    end
    println("")

    println("Side 2:")
    println("   M: ", shock["M2"])
    println("   M2n: ", shock["M2n"])
    println("   M2t: ", shock["M2t"])
    if shock["staticT"]
        println("   V2: ", shock["V2"])
        println("   V2n: ", shock["V2n"])
        println("   V2t: ", shock["V2t"])
    end
    println("")

    println("Overall Properties: ")
    println("   θ1: ", shock["thetas"][1])
    println("   θ2: ", shock["thetas"][2])
    println("   θmax: ", shock["θmax"])
    println("   δmax: ", shock["δmax"])
    println("   P2/P1: ", shock["P2P1"])
    println("   T2/T1: ", shock["T2T1"])
    println("   P02/P01: ", shock["P02P01"])
    if shock["staticT"]
        println("   T1: ", shock["T1"])
        println("   T2: ", shock["T2"])
    end
    if shock["staticP"]
        println("   P1: ", shock["P1"])
        println("   P2: ", shock["P2"])
    end
    if shock["staticT"]
        println("   a1: ", shock["a1"])
        println("   a2: ", shock["a2"])
    end
end
export printshock