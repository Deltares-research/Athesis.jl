@kernel function addToLayer!(x, dx, layerIdx)

    i, j = @index(Global, NTuple)
    @inbounds x[i,j,layerIdx] = x[i,j,layerIdx] + dx
  
end