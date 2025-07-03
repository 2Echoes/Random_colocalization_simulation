# Module for random co-localization modelisation

# volume : number of positions available for single molecules 
# abundancy : Number of single molecule in volume

module colocalization_modelisation

export self_colocalization_expectancy, self_colocalization_std, colocalization_expectancy, colocalization_std, Ncolocalization_expectancy, Ncolocalization_std


"""
The probability that voxel m is occupied with at least one molecule i.
"""
p(abundancy:: Int, volume :: Int):: Float64 = 1 - (1-1/volume)^abundancy 


"""
The probability that voxel m and l != m are both occupied with at least one molecule i
"""
q(abundancy:: Int, volume :: Int):: Float64 = (1-2*(1-1/volume)^abundancy + (1-2/volume)^abundancy)

"""
The covariance of the occupancy of two different voxel by the same distribution i.
"""
c(abundancy::Int, volume::Int)::Float64 = q(abundancy, volume) - p(abundancy,volume)^2

function self_colocalization_expectancy(abundancy::Int, volume::Int)::Float64 
    return abundancy * p(abundancy -1,volume)
end

function self_colocalization_std(abundancy::Int, volume::Int)::Float64
    p1 = p(abundancy-1, volume)

    # return sqrt(abundancy^2*p1*(1-p1))
    return sqrt(abundancy*p1*(1-p1))
end

"""
Expectancy of co-localization of distribution 1 (ie abundancy1) with distribution2 (abundancy2)
"""
function colocalization_expectancy(abundancy1::Int, abundancy2::Int, volume::Int)::Float64
    return abundancy1*p(abundancy2,volume)
end

"""
Standard deviation of co-localization of distribution 1 (ie abundancy1) with distribution2 (abundancy2)
"""
function colocalization_std(abundancy1::Int, abundancy2::Int, volume::Int)
    p2 = p(abundancy2,volume)
    c2 = c(abundancy2,volume)

    variance = abundancy1*((1-1/volume) * (p2*(1-p2) - c2) + c2)

    return sqrt(variance)
end

"""
Expectancy deviation of voxel co-occupancy with at least 1 element of all distributions (ie abundancies).
"""
function Ncolocalization_expectancy(abundancies::Vector{Int}, volume::Int)::Float64
    p_combination = prod([p(abundancy, volume) for abundancy in abundancies])
    return volume*p_combination
end

"""
Standard deviation of voxel co-occupancy with at least 1 element of all distributions (ie abundancies).
"""
function Ncolocalization_std(abundancies::Vector{Int}, volume::Int)::Float64
    p_combination = prod([p(abundancy, volume) for abundancy in abundancies])
    q_combination = prod([q(abundancy, volume) for abundancy in abundancies])
    variance = volume*p_combination*(1-p_combination) + volume*(volume-1)*(q_combination - p_combination^2)
    
    return sqrt(variance)
end

end #Module end