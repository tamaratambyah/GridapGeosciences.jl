using BenchmarkTools

x = Point(1.0,1.0)
xs = [x for i in 1:10]

bm() = ForwardMapPanel1()(x)
@benchmark bm

bm() = ∇(ForwardMapPanel1())(x)
@benchmark bm

bm() = ForwardMapPanel1()(xs)
@benchmark bm

bm() = ∇(ForwardMapPanel1())(xs)
@benchmark bm
