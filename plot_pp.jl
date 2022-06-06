using CairoMakie, AlgebraOfGraphics

function plot_pp(df,savefile)

plt = data(df) * mapping(:value) * histogram()
draw(plt)

save(savefile,plt)

end