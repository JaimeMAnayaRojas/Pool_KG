
### Plot 2 or supplementary
### Posterior probability changes LTRE
# Guppy
# nBigMatrix = 100
min_size = 4
max_size = 35
size_cen = 18.0

U= Float64(max_size)
L=Float64(min_size)
m = nBigMatrix
h = (U - L)/m
z1 =  L .+ (collect(1:m) .- 0.5) * h
vm = DataFrame(IPMs[1][2], :auto);
size_loss = Vector(mapcols(x -> LOS(x, 0), vm)[1,:])


minimum(size_loss)
z1[findall(size_loss .< 17)]


println(DataFrame(z=  z1, los= size_loss))
histogram(size_loss, bin = 20)

pF2G = plot(z1, size_loss, label = "S(z)", legend =:right, title = "a)", titleloc = :left)

vm = DataFrame(IPMs[1][3], :auto);
size_loss = Vector(mapcols(x -> LOS(x, 0), vm)[1,:])
plot!(z1, size_loss, label = "G(z,z')")


vm = DataFrame(IPMs[1][4], :auto);
size_loss = Vector(mapcols(x -> LOS(x, 0), vm)[1,:])
plot!(z1, size_loss, label = "R(z)")

ylabel!("Posterior probability (%) \n (ΔV= NK-KG> 0)")
xlabel!("Guppy size")
hspan!([90, 100],  color = :deepskyblue, alpha = 0.2, labels = false)
hspan!([0, 10],  color = :deepskyblue, alpha = 0.2, labels = false)
hspan!([95, 100],  color = :deepskyblue, alpha = 0.2, labels = false)
hspan!([0, 5],  color = :deepskyblue, alpha = 0.2, labels = false)


plot!([min_size, max_size], [50,50], c = :gray, linestyle = :dash, lab = false, linewidth = 1.5)

# Killifish

# size-speficif LTRE
# nBigMatrix = 100
min_size = 2
max_size = 110
size_cen = 18.0

U= Float64(max_size)
L=Float64(min_size)
m = nBigMatrix
h = (U - L)/m
z1 =  L .+ (collect(1:m) .- 0.5) * h
z = 0
z1 = vcat(z,z1)
vm = DataFrame(IPMs[2][2], :auto);
size_loss = Vector(mapcols(x -> LOS(x, 0), vm)[1,:])
pF2K = plot(z1, size_loss, label = false, title = "b)", titleloc = :left)

vm = DataFrame(IPMs[2][3], :auto);
size_loss = Vector(mapcols(x -> LOS(x, 0), vm)[1,:])
plot!(z1, size_loss, label =  false)


vm = DataFrame(IPMs[2][4], :auto);
size_loss = Vector(mapcols(x -> LOS(x, 0), vm)[1,:])
plot!(z1, size_loss, label = false)

ylabel!("Posterior probability (%) \n (ΔV= NK-KG> 0)")
xlabel!("Killifish size")
hspan!([90, 100],  color = :deepskyblue, alpha = 0.2, labels = false)
hspan!([0, 10],  color = :deepskyblue, alpha = 0.2, labels = false)
hspan!([95, 100],  color = :deepskyblue, alpha = 0.2, labels = false)
hspan!([0, 5],  color = :deepskyblue, alpha = 0.2, labels = false)


plot!([min_size, max_size], [50,50], c = :gray, linestyle = :dash, lab = false, linewidth = 1.5)

# l = @layout [
#     a 
# 	[grid(2,2)]
# ]


plot(pF2G, pF2K,
    layout = (2,1), size = (400, 600), titlefont = 11, foreground_color_legend = nothing, 
    margin = [2mm 3mm], xlims =(5,100)
)

