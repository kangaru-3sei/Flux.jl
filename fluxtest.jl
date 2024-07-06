using Flux

nt = (a = [2, 2], b = [2, 0], c = tanh);

g(x::NamedTuple) = sum(abs2, x.a .- x.b);

g(nt)

dg_nt = gradient(g, nt)[1]

gradient((x, y) -> sum(abs2, x.a ./ y .- x.b), nt, [1, 2])

gradient(nt, [1, 2]) do x, y
         z = x.a ./ y
         sum(abs2, z .- x.b)
       end
((a = [0.0, 0.5], b = [-0.0, -1.0], c = nothing), [-0.0, -0.25])