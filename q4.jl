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

# Solver
    model =  Model(with_optimizer(Gurobi.Optimizer))    
    
# Parameters
    
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

    # C : index set of candidate cliques 
    C = [1:amount_of_genes*amount_of_genes;]


# Variables
    @variable(model, x[1:amount_of_genes, 1:length(C)], Bin) 
    @variable(model, z[1:length(C)], Bin)
    

# Constraints 

    for i in 1:amount_of_genes
        for c in 1:length(C)
            @constraint(model, x[i,c] <= z[c])
        end 
    end 
   
    for i in 1:amount_of_genes
        @constraint(model, 1 <= sum(x[i,c] for c = 1:length(C)) <= 2)
    end

    for c in 1:length(C)
        for i in 1:amount_of_genes
            @constraint(model, sum(z[c]*x[i,c] for i = 1:amount_of_genes) >= 2*z[c])
        end 
    end 


    for i in 1:amount_of_genes
        for j in 1:amount_of_genes
            for c in 1:length(C)
                if i != j

                    if adjacency((genes[i]),(genes[j])) != 0 
                        not_edge = 0
                    else 
                        not_edge = 1
                    end 

                    @constraint(model, not_edge * (x[i,c] + x[j,c]) <= 1)

                end
            end
        end
    end


# Objective function
    @objective(model, Min, sum(z[i] for i = 1:length(C)))

# Invoke solver
    optimize!(model)

# Display results
    println("Objective value: ---> ", JuMP.objective_value(model))
