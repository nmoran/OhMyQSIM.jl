
export reshape_tensor

function reshape_tensor(tensor, indices)
    index_order = Array{Int, 1}()
    dims = Array{Int, 1}()
    for x in indices
        dim = 1
        for y in x
            push!(index_order, y)
            dim *= size(tensor)[y]
        end
        push!(dims, dim)
    end
    tensor = permutedims(tensor, index_order)
    reshape(tensor, Tuple(dims))
end
