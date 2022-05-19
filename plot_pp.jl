using CairoMakie, AlgebraOfGraphics

function plot_pp(df,savefile)

plt = data(df) * mapping(:valued) * histogram()
draw(plt)

save(savefile,plt)

end