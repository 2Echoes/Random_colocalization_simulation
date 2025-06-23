module coordinates_generator

export generate_coordinates



"""
Genereate random coordinates within boundaries given in `shape` parameters. Output is shaped (`simulation_number`, `coordinates_number` , dim(`shape`)).
"""
function generate_coordinates(shape::Tuple, coordinates_number::Int, simulation_number=1::Int)::Array{Int,3}
    return [rand(1:dim) for simulation in 1:simulation_number, coordinate in 1:coordinates_number, dim in shape]
end

end