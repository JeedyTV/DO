using JuMP
using Gurobi

# Declare the solver 
    model =  Model(with_optimizer(Gurobi.Optimizer))   

# Define the parameters
    file =  open("Hsa/BGSE1456.txt", "r")
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
      
# Create the variables 
    @variable(model, x[1:amount_of_genes], Bin)

# Define the constraints 
    for i in 1:amount_of_genes
        for j in 1:amount_of_genes
            if i != j
                if (genes[i], genes[j]) ∉ coexpression && (genes[j], genes[i]) ∉ coexpression
                    @constraint(model,(x[i] + x[j]) <= 1)
                end 
            end
        end
    end

# Define the objective function
    @objective(model, Max, sum(x[i] for i = 1:amount_of_genes))

# Invoke the solver and display the result 
    optimize!(model)

# Display the results
    println("Objective value : ---> ", JuMP.objective_value(model))
    
