# 確率密度の比較のための plot
# バックエンド: pyplot

using Distributions
using Plots
pyplot()

x = collect(-100:0.1:100)
y = pdf.(Normal(3,10),x)
z = rand(Normal(3,10), 100)

histogram(z, nbins = 20, normed=true)
plot!(x,y)


# 各 chain ごとの (変数の値,iteration number)のプロット
# バックエンド: pyplot
