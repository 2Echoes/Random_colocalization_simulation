module colocalization
# Module to quantify colocalization between arrays of coordinates

using Distances
using LinearAlgebra

export colocalization_test, count_colocalization, compute_colocalization_across_simulations, compute_colocalization_count_across_simulations




"""
Return the distances from query points to anchor_points. 
"""
function _find_closest_coordinates(anchor_points::AbstractArray{<:Int,2}, query_points::AbstractArray{<:Int,2})
    distance_matrix = pairwise(Euclidean(), anchor_points', query_points')
    distance_to_closest_neighbor = findmin(eachrow(distance_matrix))[1]
    return distance_to_closest_neighbor
end

"""
Perform a co-localization test between two array of coordinates. For each coordinates in `query_points` return 1 if a coordinate closer than `distance_threshold` found
in `anchor_points`.
"""
function colocalization_test(
    anchor_points::AbstractArray{<:Int,2}, 
    query_points::AbstractArray{<:Int,2}, 
    distance_threshold::Int
    )::Vector{Bool}
    
    distances_to_neighbor = _find_closest_coordinates(anchor_points, query_points)

    return distances_to_neighbor .<= distance_threshold

end

function count_colocalization(
    anchor_points::AbstractArray{<:Int,2}, 
    query_points::AbstractArray{<:Int,2}, 
    distance_threshold::Int
    ) :: Int

    return sum(colocalization_test(anchor_points, query_points, distance_threshold))
end

"""
Compute colocalization_test across a 3D array where 3rd dim is representing simulations.
"""
function compute_colocalization_across_simulations(
    anchor_points::Array{Int,3},
    query_points::Array{Int,3},
    distance_threshold::Int
    ) :: Array{Bool,2}

    simulation_number, distribution_number,_ = size(query_points)
    result = Array{Bool,2}(undef, simulation_number, distribution_number)

    for simulation in 1:simulation_number
        result[simulation,:] = colocalization_test(
            view(anchor_points, simulation, :, :),
            view(query_points, simulation, :, :),
            distance_threshold
        )
    end

    return result
end

function compute_colocalization_count_across_simulations(
    anchor_points::Array{Int,3},
    query_points::Array{Int,3},
    distance_threshold::Int
    )

    result = map(
        count_colocalization, 
        eachslice(anchor_points,dims=1), 
        eachslice(query_points,dims=1), 
        fill(distance_threshold, size(anchor_points,1))
    )
    return result
end


end #end module