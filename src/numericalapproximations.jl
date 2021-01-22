
function DivergingFlowAreaMethod(Lvec, Avec, Pb, Patm, Tstar, gamma;f=0.04, n=100, shape="square", R=287, μ=1.81e-5)
    ### Fit Length and Area distribution
    Lfit = Akima(Lvec, Avec)
    if typeof(f)==Array{Float64,2}
        ffit = Akima(f[:,1], f[:,2])
    else
        ffit(Re) = f
    end

    ### Create nodes (Adistro and Ldistro)
    Ldistro = collect(range(Lvec[1], Lvec[end], length=n))
    Adistro = Lfit.(Ldistro)

    ### Block 1 Isentropic Prepwork
    Mnvec = zeros(n)
    Pnvec = zeros(n)
    PnP0vec = zeros(n)

    Pe=Pb
    PePstar = Pe/Patm
    Me = PMfun(gamma, PePstar)
    AeAstar = Afun(gamma, Me)
    Ae = Adistro[1]
    An = Ae

    Mn = Me
    Pn = Pe
    AnAstar = AeAstar

    TnTstar = Tfun(gamma, Mn)
    Tn = TnTstar*Tstar
    an = afun(gamma, R, Tn)
    Vn = Mn*an
    ρn = Pn/(R*Tn)
    dn = sqrt(An)
    Ren = ρn*Vn*dn/μ
    fn = ffit(Ren)

    Mnvec[1] = Mn
    Pnvec[1] = Pn
    PnP0vec[1] = PePstar
    

    for i=1:n-1
        ### Block 2 - Fanno Flow
        Ln = Ldistro[i]-Ldistro[i+1] #Should be the same time to time because I linearly spaced Ldistro. 
        # println("Ln: ", Ln)
        An = Adistro[i]
        if shape=="square"
            dn = sqrt(An)
        else shape=="circle"
            dn = sqrt(4*An/pi)
        end
        # println("fn: ", fn)
        fld = fn*Ln/dn
        fld2 = fldfun(Mn, gamma)
        fld1 = fld + fld2
        M1 = fldMfun(fld1, gamma)
        ### Block 3 - Isentropic flow
        A1Astar = Afun(gamma, M1)
        P1Pstar = Pfun(gamma, M1)
        Astar = (An/A1Astar) #*0.999
        An1 = Adistro[i+1] #*1.001
        An1Astar = An1/Astar
        # println(An1Astar)
        Mn1 = 0
        try
            Mn1 = AMfun(gamma, An1Astar) 
        catch
            break
        end
        # println(Mn1)
        Pn1Pstar = Pfun(gamma, Mn1)
        Pstar = Pn/P1Pstar
        Pn1 = Pn1Pstar*Pstar
        # println("Pn: ", Pn)

        Mnvec[i+1] = Mn1
        PnP0vec[i+1] = Pn1Pstar
        Pnvec[i+1] = Pn1 #Need to calculate Pstar

        # #Testing, this is what it should be just being isentropic flow. 
        # Atstar = Ae/AeAstar
        # An1Atstar = An1/Atstar
        # Mn1o = AMfun(gamma, An1Atstar)
        # Pn1oPstar = Pfun(gamma, Mn1o)
        
        # println(Pn)
        # println("Run $i") #I was overwriting Mn before printing ... Ooops
        # println("Mn: ", Mn)
        # println("M1: ", M1)
        # println("Mn1: ", Mn1)
        # println("Mn1o: ", Mn1o)
        # println("fanno speed diff: ", M1-Mn)
        # println("tot speed diff: ", Mn1-Mn)
        # println("isen speed diff: ", Mn1o-Mn)
        # println("comp speed diff: ", Mn1-Mn1o)
        # println("")
        # println("Pn1Pstar: ", Pn1Pstar)
        # println("Pn1oPstar: ", Pn1oPstar)
        # println("press diff: ", Pn1Pstar- Pn1oPstar)
        # println("")
        # println("An: ", An)
        # println("A1: ", A1)
        # println("An1: ", An1)
        # println("A1Astar: ", A1Astar)
        # println("")

        ### Move to next node
        Mn = Mn1
        AnAstar = An1Astar
        Pn = Pn1

        TnTstar = Tfun(gamma, Mn)
        Tn = TnTstar*Tstar
        an = afun(gamma, R, Tn)
        Vn = Mn*an
        ρn = Pn/(R*Tn)
        dn = sqrt(An1)
        Ren = ρn*Vn*dn/μ
        fn = ffit(Ren)
    end
    outs = Dict()
    outs["Mnvec"] = Mnvec
    outs["Pnvec"] = Pnvec
    outs["PnP0vec"] = PnP0vec
    outs["Ldistro"] = Ldistro
    return outs
end
export DivergingFlowAreaMethod


function DivergingFlowForward(Lvec, Avec, Pb, Patm, Tstar, gamma;f=0.04, n=100, shape="square", R=287, μ=1.81e-5, subsonic=true)
    ### Fit Length and Area distribution
    Lfit = Akima(Lvec, Avec)
    if typeof(f)==Array{Float64,2}
        ffit = Akima(f[:,1], f[:,2])
    else
        ffit(Re) = f
    end

    ### Create nodes (Adistro and Ldistro)
    Ldistro = collect(range(Lvec[1], Lvec[end], length=n))
    Adistro = Lfit.(Ldistro)

    ### Block 1 Isentropic Prepwork
    Mnvec = zeros(n)
    Pnvec = zeros(n)
    PnP0vec = zeros(n)

    Pe=Pb
    PePstar = Pe/Patm
    Me = PMfun(gamma, PePstar; subsonic=subsonic)
    AeAstar = Afun(gamma, Me)
    Ae = Adistro[1]
    An = Ae

    Mn = Me
    Pn = Pe
    AnAstar = AeAstar

    TnTstar = Tfun(gamma, Mn)
    Tn = TnTstar*Tstar
    an = afun(gamma, R, Tn)
    Vn = Mn*an
    ρn = Pn/(R*Tn)
    dn = sqrt(An)
    Ren = ρn*Vn*dn/μ
    fn = ffit(Ren)

    Mnvec[1] = Mn
    Pnvec[1] = Pn
    PnP0vec[1] = PePstar
    

    for i=1:n-1
        ### Block 2 - Fanno Flow
        Ln = Ldistro[i+1]-Ldistro[i] #Should be the same time to time because I linearly spaced Ldistro. 
        # println("Ln: ", Ln)
        An = Adistro[i]
        if shape=="square"
            dn = sqrt(An)
        else shape=="circle"
            dn = sqrt(4*An/pi)
        end
        fld = fn*Ln/dn
        fld1 = fldfun(Mn, gamma)
        fld2 = fld1 - fld
        # println("Mn: ", Mn)
        # println("fld: ", fld)
        # println("fld1: ", fld1)
        # println("fld2: ", fld2)
        try
            M2 = fldMfun(fld2, gamma; subsonic=false)
        catch
            break
        end
        # println(M2)
        ### Block 3 - Isentropic flow
        A2Astar = Afun(gamma, M2)
        P2Pstar = Pfun(gamma, M2)
        Astar = (An/A2Astar) #*0.999
        An1 = Adistro[i+1] #*1.001
        An1Astar = An1/Astar
        Mn1 = AMfun(gamma, An1Astar; subsonic=subsonic) 
        Pn1Pstar = Pfun(gamma, Mn1)
        Pstar = Pn/P2Pstar
        Pn1 = Pn1Pstar*Pstar
        # println("Pn: ", Pn)

        Mnvec[i+1] = Mn1
        PnP0vec[i+1] = Pn1Pstar
        Pnvec[i+1] = Pn1 #Need to calculate Pstar

        # #Testing, this is what it should be just being isentropic flow. 
        Atstar = Ae/AeAstar
        An1Atstar = An1/Atstar
        Mn1o = AMfun(gamma, An1Atstar; subsonic=subsonic)
        Pn1oPstar = Pfun(gamma, Mn1o)
        
        # println(Pn)
        # println("Run $i") #I was overwriting Mn before printing ... Ooops
        # println("Mn: ", Mn)
        # println("M1: ", M1)
        # println("Mn1: ", Mn1)
        # println("Mn1o: ", Mn1o)
        # println("fanno speed diff: ", M1-Mn)
        # println("tot speed diff: ", Mn1-Mn)
        # println("isen speed diff: ", Mn1o-Mn)
        # println("comp speed diff: ", Mn1-Mn1o)
        # println("")
        # println("Pn1Pstar: ", Pn1Pstar)
        # println("Pn1oPstar: ", Pn1oPstar)
        # println("press diff: ", Pn1Pstar- Pn1oPstar)
        # println("")
        # println("An: ", An)
        # println("A1: ", A1)
        # println("An1: ", An1)
        # println("A1Astar: ", A1Astar)
        # println("")

        ### Move to next node
        Mn = Mn1
        AnAstar = An1Astar
        Pn = Pn1

        TnTstar = Tfun(gamma, Mn)
        Tn = TnTstar*Tstar
        an = afun(gamma, R, Tn)
        Vn = Mn*an
        ρn = Pn/(R*Tn)
        dn = sqrt(An1)
        Ren = ρn*Vn*dn/μ
        fn = ffit(Ren)
    end
    outs = Dict()
    outs["Mnvec"] = Mnvec
    outs["Pnvec"] = Pnvec
    outs["PnP0vec"] = PnP0vec
    outs["Ldistro"] = Ldistro
    return outs
end
export DivergingFlowForward




# function ComplexDiff(fun, x)
# # I haven't checked this yet, turns out I didn't need it for the project because I had an analytical function for my Area... so I could get the analytical derivative. 
#     Alpha = 1e-8
#     Alphai = Alpha*im
#     yprime = zeros(length(x))
#     for i = 1:length(x)
#         xprime = x
#         xprime[i] = x[i]+Alphai
#         yprime[i] = imag(fun.(xprime))/Alpha
#     end
#     return yprime
# end
# export ComplexDiff
"""
    RK4(xi, yi, Δx, fun)

Runge-Kutta 4th ODE solver. 

### Inputs
- xi::Float64 - Current x location
- yi::Float64 - Current y location
- Δx::Float64 - x distance from current node to next node
- fun::Function Handle of derivative function to be integrated. 
"""
function RK4(xi, yi, Δx, fun)
    k1 = fun(xi, yi)
    k2 = fun(xi+(Δx/2), yi+(k1*Δx/2))
    k3 = fun(xi+(Δx/2), yi+(k2*Δx/2))
    k4 = fun(xi+Δx, yi+(k3*Δx))
    yi1 = yi + Δx*(k1 + 2*k2 + 2*k3 + k4)/6
    return yi1
end
export RK4

function nearestto(x, xp)
    residuals = abs.(x .-xp)
    min, idx = findmin(residuals)
    return x[idx], idx
end
export nearestto

# function findzero(fun, rng)
#     vector = collect(rng[1]:.00001:rng[2])
#     result = zeros(length(vector))
#     for i=1:length(vector)
#         result[i] = fun(vector[i])
#     end
#     zero, idx = nearestto(result, 0)
#     return vector[idx]
# end

function errfun(x, xt)
    return 100*(x-xt)/xt
end
export errfun