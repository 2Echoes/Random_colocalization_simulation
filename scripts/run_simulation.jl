# Colocalization simulation.jl

# author : Floric Slimani(floric.slimani@igh.cnrs.fr)

#Import
include("../src/coordinates_generation.jl")
include("../src/colocalization.jl")

using ProgressMeter
using .coordinates_generator
using .colocalization
using DataFrames
using JLD2
using CSV

#Simulation Parameters
simulation_number = 100000
distance_threshold = 0
shape = (1, 158,158)
distributions_abundacy = [10,50,100,200,250,400,500]

"""
Performs a random co-localization simulation and saves results at the output_path parameters with name 'colocalization_counts_table' which will be an Int array of size (simulation_number,distributions_number,distributions_number)
Self co-localization are not corrected.
"""
function main(
    output_path =               length(ARGS) >=1 && ARGS[1] !='_' ? ARGS[1] : pwd(),
    simulation_number =         length(ARGS) >=2 && ARGS[2] !='_' ? ARGS[2] : simulation_number,
    distance_threshold =        length(ARGS) >=3 && ARGS[3] !='_' ? ARGS[3] : distance_threshold,
    shape =                     length(ARGS) >=4 && ARGS[4] !='_' ? ARGS[4] : shape,
    distributions_abundacy =    length(ARGS) >=5 && ARGS[5] !='_' ? ARGS[5] : distributions_abundacy
    )
    
    #Init variables
    distributions_coordinates = Dict()
    distributions_number=length(distributions_abundacy)
    distributions_abundacy_dict = Dict(zip(1:distributions_number, distributions_abundacy))
    
    #Sum up of user parameters
    Sumup_df = DataFrame(
        distribution_id = collect(keys(distributions_abundacy_dict)), 
        abudancy= collect(values(distributions_abundacy_dict)), 
        )
        CSV.write(output_path*"/colocalization_simulation_parameters.csv", Sumup_df)
        
    println("Starting $simulation_number simulations for $distributions_number distributions :")
    
    #Generate coordinates
    print("Generating random coordinates :")
    @time for (key,value) in distributions_abundacy_dict
        coordinates = generate_coordinates(shape, value, simulation_number)
        distributions_coordinates[key] = coordinates
    end

    #Test pair colocalization
    println("Counting pair co-localization :")
    colocalization_counts = Array{Int}(undef,simulation_number,distributions_number,distributions_number)
    @showprogress for line in 1:distributions_number
        for col in 1:distributions_number

            colocalization_counts[:,line,col] = compute_colocalization_count_across_simulations(
                                                distributions_coordinates[col],
                                                distributions_coordinates[line],
                                                distance_threshold
                                                )
        end
    end

    #Building truth table for all co-localization
    println("Testing all spots for all distribututions :")
    total_number_spots = sum(distributions_abundacy)
    colocalization_truth_table = Array{Bool,3}(undef,simulation_number,total_number_spots,distributions_number)
    colocalization_truth_table_index = vcat(
        (fill(idx, distributions_abundacy[idx]) for idx in 1:distributions_number)...
     ) # index for colocalization_counts_table whith consistant order and value representing which distribution the line of truth table belongs to.
    
    
    ## Initializing index for loop 
    spot_idx_inf = 1
    spot_idx_sup = nothing
    @showprogress for line in 1:distributions_number
        if line != 1 
            spot_idx_inf = spot_idx_sup + 1
        end
        spot_idx_sup = spot_idx_inf + distributions_abundacy[line] -1
        for col in 1:distributions_number
            colocalization_truth_table[:,spot_idx_inf : spot_idx_sup,col] = compute_colocalization_across_simulations(
            distributions_coordinates[col],
            distributions_coordinates[line],
            distance_threshold
            )
        end
    end


    println("Saving results at $output_path")
    @save output_path*"/colocalization_pairwise_counts.jld2" colocalization_counts
    @save output_path*"/colocalization_truth_table.jld2" colocalization_truth_table colocalization_truth_table_index shape

end #end main

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end