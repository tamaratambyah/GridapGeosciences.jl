
function plot_error(ns,errs;
  leginf = fill(false,Int(length(errs)/length(ns))),
  ls=[:solid, :dash, :dot, :dashdot, :dashdotdot],
  colors = palette(:tab10),
  markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross],
  ms=[6,6,8,6,6,8,8] )

  nsims = Int(length(errs)/length(ns))
  for i in 1:nsims
    idx1 = 1 + (i-1)*length(ns)
    idx2 = (i)*length(ns)
    plot!(ns,
          errs[idx1:idx2],
          lw=3,
          markersize=ms[i],
          c=colors[i],ls=ls[i], markershape=markers[i],
          label=leginf[i])
  end

end
