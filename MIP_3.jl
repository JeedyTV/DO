using JuMP
using Gurobi

# function which represents the adjacency matrix
function adjacency(gene1, gene2)
    if (gene1, gene2) in coexpression || (gene2, gene1) in coexpression
        return 1
    else 
        return 0
    end
end

# Declare the solver 
    model =  Model(with_optimizer(Gurobi.Optimizer))   

# Define the parameters
    
    #not functional because too long
    file =  open("Hsa/BGSE1456.txt", "r")
    
    #functional file, simpler graph
    #file =  open("simpleGraph.txt", "r")

    data = readlines(file);

    # Create a list with all the gene labels
    genes = []
    coexpression = []
    for line in data
        gene1, gene2, c = split(line, '\t')
        append!(coexpression, [(gene1, gene2)])
        append!(genes, [gene1, gene2])     
    end
    genes = unique(genes) 
    amount_of_genes = length(genes)
    # julia set has O(1) access
    coexpression = Set(coexpression)

    #density = 1
    density = 10/15

# Create the variables 
    @variable(model, x[1:amount_of_genes], Bin)
    @variable(model, z[1:amount_of_genes , 1:amount_of_genes], Int)

# Define the constraints 
    @constraint(model, sum( (sum( (z[i,j] * ( density - adjacency(genes[i],genes[j]) )) for j = i+1:amount_of_genes)) for i in 1:amount_of_genes-1 ) <= 0)

    for i in 1:amount_of_genes-1 
        for j in i:amount_of_genes
            @constraint(model, z[i,j] <= x[i])
            @constraint(model, z[i,j] <= x[j])
            @constraint(model, z[i,j] + 1 >= x[i] + x[j])
            @constraint(model, z[i,j] >= 0)
        end
    end

# Define the objective function
    @objective(model, Max, sum(x[i] for i = 1:amount_of_genes))

# Invoke the solver and display the result 
    optimize!(model)

# Display the results
    println("Objective value : ---> ", JuMP.objective_value(model))
    