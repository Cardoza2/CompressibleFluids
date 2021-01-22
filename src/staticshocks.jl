export normshockTfun, normshockPfun, normshockrhofun, normshockMfun, normshockPnfun, normshockPnMfun, normshockP0fun, normshock, eq426, findshockarea

function normshockTfun(gamma, M)
    topleft = 1 + (((gamma-1)/2)*(M^2))
    topright = ((2*gamma)/(gamma-1))*(M^2) -1
    bottom = (((gamma+1)^2)/(2*(gamma-1)))*M^2
    return topleft*topright/bottom
end

function normshockPfun(gamma, M)
    top1 = 2*gamma*(M^2)
    bot1 = gamma + 1
    top2 = gamma-1
    bot2 = gamma+1
    return (top1/bot1)-(top2/bot2)
end

function normshockrhofun(gamma, M)
    top = (gamma+1)*(M^2)
    bot = (gamma-1)*(M^2) + 2
    return top/bot
end

function normshockMfun(gamma, M)
    top = (M^2) + (2/(gamma-1))
    bot = (2*gamma*(M^2))/(gamma-1) - 1
    return sqrt(top/bot)
end

function normshockPnfun(gamma, M)
    topleft = (M^2)*(gamma+1)/2
    botleft = (1 + (M^2)*(gamma-1)/2)
    botright = (M^2)*2*gamma/(gamma+1) - (gamma-1)/(gamma+1)
    return ((topleft/botleft)^(gamma/(gamma-1)))*((1/botright)^(1/(gamma-1)))
end

function normshockPnMfun(Gamma, P02P01)
    rng = (1, 30)
    fun(m; gamma=Gamma, p02p01=P02P01) = normshockPnfun(gamma, m) - p02p01
    mach = find_zero(fun, rng, Bisection())
    return mach
end

function normshockPMfun(Gamma, P2P1)
    rng = (1, 30)
    fun(m; gamma=Gamma, p2p1=P2P1) = normshockPfun(gamma, m) - p2p1
    mach = find_zero(fun, rng, Bisection()) #These find_zero functions could be replaced with a real root finder. I think fairly easily too. 
    return mach
end
export normshockPMfun

function normshockP0fun(gamma, M)
    temp = (1 + (M^2)*(gamma-1)/2)^(gamma/(gamma-1))
    return 1/temp
end

function normshock(gamma, M)
    shock = Dict()
    shock["M2"] = normshockMfun(gamma, M)
    shock["T2T1"] = normshockTfun(gamma, M)
    shock["P2P1"] = normshockPfun(gamma, M)
    shock["rho2rho1"] = normshockrhofun(gamma, M)
    shock["P02P01"] = normshockPnfun(gamma, M)
    shock["P1P01"] = normshockP0fun(gamma, M)
    shock["A1starA2star"] = shock["P02P01"]
    shock["P2P2star"] = Pfun(gamma, shock["M2"])
    return shock
end

function eq426(PbP01, AeAt, gamma)
    #Use equation 4.26 to calculate Me
    t0 = -1/(gamma-1)
    t1 = (1/(gamma-1))^2
    t2 = 2/(gamma-1)
    t3 = (2/(gamma+1))^((gamma+1)/(gamma-1))
    t4 = (1/PbP01)^2
    t5 = (1/AeAt)^2
    Me2 = t0 + sqrt(t1 + t2*t3*t4*t5)
    return sqrt(Me2)
end

function findshockarea(AeAt, PbP01, gamma)
    Me = eq426(PbP01, AeAt, gamma) #Find the exit mach using equation 4.26
    Pe_P02 = Pfun(gamma, Me)  #Use that mach number to find the Pe pressure ratio
    P02P01 = PbP01/Pe_P02 #Find the shock stagnation pressure ratio
    M1 = normshockPnMfun(gamma, P02P01) #Use that relationship to find M1
    # println("M1: ", M1)
    A1A1star = Afun(gamma, M1) #Use M1 to find A1A1star
    return A1A1star
end