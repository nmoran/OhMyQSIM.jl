
export reshape_tensor, tensor_identity

"""
    function reshape_tensor(tensor, indices)

Permutes and reshapes tensor to move and merge indices
"""
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

function tensor_identity(m, n, d)
    id = zeros(m, n, d, d)
    id_small = Diagonal(ones(d))
    for i in 1:m
        for j in 1:n
            if i == j
                id[i, j, :, :] = id_small
            end
        end
    end
    id
end
