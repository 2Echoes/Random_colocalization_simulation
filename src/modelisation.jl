# Module for random co-localization modelisation

# volume : number of positions available for single molecules 
# abundancy : Number of single molecule in volume

module colocalization_modelisation

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
c(abundacy::Int, volume::Int)::Float64 = q(abundacy, volume) - p(abundacy,volume)^2

export self_colocalization_expectancy
function self_colocalization_expectancy(abundancy::Int, volume::Int)::Float64 
    return abundancy - volume * p(abundancy,volume)
end

export self_colocalization_std
function self_colocalization_std(abundancy::Int, volume::Int)::Float64
    p1 = p(abundancy, volume)
    c1 = c(abundancy, volume)

    return sqrt(volume*p1(1-p1) + volume*(volume-1)*c1)
end

export colocalization_expectancy
"""
Expectancy of co-localization of distribution 1 (ie abundancy1) with distribution2 (abundacy2)
"""
function colocalization_expectancy(abundacy1::Int, abundacy2::Int, volume::Int)::Float64
    return (abundacy1/volume)*p(abundacy2,volume)
end

export colocalization_std
"""
Standard deviation of co-localization of distribution 1 (ie abundancy1) with distribution2 (abundacy2)
"""
function colocalization_std(abundacy1::Int, abundacy2::Int, volume::Int)
    p2 = p(abundacy2,volume)
    c2 = c(abundacy2,volume)

    variance = abundacy1*((1-1/volume) * (p2*(1-p2) - c2) + c2)

    return sqrt(variance)
end

export Ncolocalization_expectancy
"""
Expectancy deviation of voxel co-occupancy with at least 1 element of all distributions (ie abundancies).
"""
function Ncolocalization_expectancy(abundancies::Vector{Int}, volume::Int)::Float64
    p_combination = prod([p(abundacy, volume) for abundancy in abundancies])
    return v*p_combination
end

export Ncolocalization_std
"""
Standard deviation of voxel co-occupancy with at least 1 element of all distributions (ie abundancies).
"""
function Ncolocalization_std(abundancies::Vector{Int}, volume::Int)::Float64
    p_combination = prod([p(abundacy, volume) for abundancy in abundancies])
    q_combination = prod([q(abundacy, volume) for abundancy in abundancies])
    variance = v*p_combination*(1-p_combination) + v(v-1)*(q_combination - p_combination^2)
    
    return sqrt(variance)
end

end #Module end