function norm(A)
    temp = 0
    for i=1:length(A)
        temp += A[i]^2
    end
    return sqrt(temp)
end
export norm


"""
#### asatriangle(a, C, b)

### Inputs
- a - the angle adjacent to leg A
- C - the length of the bottom of the triangle
- b - the angle adjacent to leg B

### Notes
I made a booboo, and define angle a as the angle adjacent to A, which is not what is typically done in math, but I don't want to go back and fix that right now. 
Assuming that angles a = b. -> need to adjust this so that I can have different angled triangles, but I doubt it'll be on the test. 
"""
function asatriangle(a, C, b; α=0, showplot=false)
    triangle = Dict()
    triangle["a"] = a
    triangle["b"] = b
    triangle["c"] = 180-a-b
    triangle["C"] = C
    triangle["H"] = C*tand(a)/2
    Rmat = [cosd(-α) -sind(-α); sind(-α) cosd(-α)]
    triangle["Avec"] = Rmat*[C/2; triangle["H"]]
    triangle["Bvec"] = Rmat*[C/2; -triangle["H"]]
    triangle["Cvec"] = Rmat*[-C; 0]
    triangle["A"] = norm(triangle["Avec"])
    triangle["B"] = norm(triangle["Bvec"])
    triangle["Δa"] = 180-triangle["c"]
    triangle["Δb"] = 180-triangle["c"]
    triangle["astar"] = 180-a
    triangle["bstar"] = 180-b
    triangle["Anorm"] = Rmat*([triangle["H"]; -C/2]/norm([triangle["H"]; -C/2]))
    triangle["Bnorm"] = Rmat*([-triangle["H"]; -C/2]/norm([-triangle["H"]; -C/2]))
    triangle["Cnorm"] = Rmat*[0; 1]
    triangle["α"] = α

    #Prepare the plot
    p1 = [0; 0]
    p2 = p1 + [C/2; triangle["H"]]
    p3 = p2 + [C/2; -triangle["H"]]
    p4 = p3 + [-C; 0]
    p1 = Rmat*p1
    p2 = Rmat*p2
    p3 = Rmat*p3
    p4 = Rmat*p4
    xy = hcat(p1, p2, p3, p4)
    triangle["indicies"] = xy[:,1:3]
    Ap = p1 + Rmat*([C/2; triangle["H"]]./2) + [0; p2[2]*.1]
    Bp = p2 + Rmat*([C/2; -triangle["H"]]./2) + [0; p2[2]*.1]
    Cp = p3 + Rmat*([-C; 0]./2) + [0; p2[2]*.1]

    # println(xy)
    triplt = plot(xy[1,:], xy[2,:], linecolor=:black, linewidth=4, leg=false, lab="Triangle")
    annotate!([(Ap[1], Ap[2], Plots.text("A", 10, :black, :center)), (Bp[1], Bp[2], Plots.text("B", 10, :black, :center)), (Cp[1], Cp[2], Plots.text("C", 10, :black, :center))])
    triangle["plot"] = triplt

    if showplot
        display(triplt)
        println("Plots look distorted because of scaling.")
    end
    return triangle
end
export asatriangle

function shtriangle(C, H; α=0, showplot=false)
    triangle = Dict()
    triangle["C"] = C
    triangle["H"] = H
    D = C/2
    triangle["a"] = atand(H/D)
    triangle["b"] = triangle["a"]
    triangle["c"] = 180 - triangle["a"] - triangle["b"]
    Rmat = [cosd(-α) -sind(-α); sind(-α) cosd(-α)]
    triangle["Avec"] = Rmat*[D; H]
    triangle["Bvec"] = Rmat*[D; -H]
    triangle["Cvec"] = Rmat*[-C; 0]
    triangle["A"] = norm(triangle["Avec"])
    triangle["B"] = norm(triangle["Bvec"])
    triangle["Anorm"] = Rmat*([H; -D]./norm(triangle["Avec"]))
    triangle["Bnorm"] = Rmat*([-H; -D]./norm(triangle["Avec"]))
    triangle["Cnorm"] = [0; 1]
    triangle["Δa"] = 180-triangle["c"]
    triangle["Δb"] = 180-triangle["c"]
    triangle["astar"] = 180-triangle["a"]
    triangle["bstar"] = 180-triangle["b"]
    triangle["α"] = α

    #Prepare the plot
    p1 = [0; 0]
    p2 = p1 + [D; H]
    p3 = p2 + [D; -H]
    p4 = p3 + [-C; 0]
    p1 = Rmat*p1
    p2 = Rmat*p2
    p3 = Rmat*p3
    p4 = Rmat*p4
    xy = hcat(p1, p2, p3, p4)
    triangle["indicies"] = xy[:,1:3]
    Ap = p1 + Rmat*([D; H]./2) + [0; p2[2]*.1]
    Bp = p2 + Rmat*([D; -H]./2) + [0; p2[2]*.1]
    Cp = p3 + Rmat*([-C; 0]./2) + [0; p2[2]*.1]

    triplt = plot(xy[1,:], xy[2,:], linecolor=:black, linewidth=4, leg=false, lab="Triangle")
    annotate!([(Ap[1], Ap[2], Plots.text("A", 10, :black, :center)), (Bp[1], Bp[2], Plots.text("B", 10, :black, :center)), (Cp[1], Cp[2], Plots.text("C", 10, :black, :center))])
    triangle["plot"] = triplt

    if showplot
        display(triplt)
        println("Plots look distorted because of scaling.")
    end
    return triangle
end
export shtriangle

