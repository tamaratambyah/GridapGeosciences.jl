function writevtk_ambient(ambient_model,
  panel_cfs::AbstractArray{<:CellField},labels::AbstractArray{<:String},
  panel_ids::AbstractArray{Int})

  Ω_sphere = Triangulation(ambient_model)
  ambient_cfs = map(x -> ambient_cellfield(x,Ω_sphere,panel_ids),panel_cfs)
  cellfields = map((x,y) -> x=>y, labels,ambient_cfs)

  writevtk(Ω_sphere,dir*"/ambient_model",
            cellfields=cellfields,
            append=false)
end


function writevtk_panel(panel_model,
  panel_cfs::AbstractArray{<:CellField},labels::AbstractArray{<:String},
  panel_ids::AbstractArray{Int})

  cellfields = map((x,y) -> x=>y, labels,panel_cfs)

  for p in collect(1:6)
    mask = panel_ids.==p
    Ωp = Triangulation(panel_model,mask)
    writevtk(Ωp,dir*"/panel$(p)_model",
            cellfields=cellfields,
            append=false)
  end


end
