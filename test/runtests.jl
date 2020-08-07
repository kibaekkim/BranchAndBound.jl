using Test
using BranchAndBound
using JuMP, GLPK
using MathOptInterface

const BB = BranchAndBound
const MOI = MathOptInterface

model = Model(GLPK.Optimizer)
@variable(model, 0 <= x <= 1)
@objective(model, Min, x)

tree = BB.AbstractTree(model, GLPK.Optimizer)
node = BB.AbstractNode()

@testset "Initial Tree" begin
    @test tree.node_counter == 0
    @test length(tree.nodes) == 0
    @test length(tree.processed) == 0
end

@testset "Initial Node" begin
    @test node.id == -1
    @test isnothing(node.parent)
    @test node.depth == 0
    @test isnothing(node.model)
    @test isnothing(node.reference_map)
    @test length(node.variable_lower_bound_changes) == 0
    @test length(node.variable_upper_bound_changes) == 0
    @test length(node.constraint_changes) == 0
    @test node.solution_status == MOI.OPTIMIZE_NOT_CALLED
    @test isinf(node.primal_bound)
    @test isinf(-node.dual_bound)
end

@testset "Tree with root" begin
    BB.push!(tree, node)
    @test tree.node_counter == 1
    @test length(tree.nodes) == 1
    @test length(tree.processed) == 0
    @test tree.nodes[1] == node
end

current_node = BB.next_node(tree)

@testset "Root node with model" begin
    @test !isnothing(current_node.model)
    @test JuMP.lower_bound(current_node.reference_map[x]) == 0.
    @test JuMP.upper_bound(current_node.reference_map[x]) == 1.
end

@testset "Bounding root node" begin
    BB.bound!(current_node)
    @test current_node.solution_status == MOI.OPTIMAL
    @test current_node.dual_bound == 0.0
    @test BB.node_solution(current_node, x) == 0.0
end

@testset "Marked root node as processed" begin
    BB.processed!(tree, current_node)
    @test isempty(tree.nodes)
    @test length(tree.processed) == 1
    @test tree.processed[1] == current_node
end

@testset "Branching at root" begin
    child1 = BB.AbstractNode(current_node)
    @test child1.parent == current_node
    @test child1.depth == 1
    @test child1.dual_bound == current_node.dual_bound
    child1.variable_lower_bound_changes[x] = 0.5
    child1.dual_bound = 0.5

    child2 = BB.AbstractNode(current_node)
    @test child2.parent == current_node
    @test child2.depth == 1
    @test child2.dual_bound == current_node.dual_bound
    child2.variable_upper_bound_changes[x] = 0.5
    child2.dual_bound = 0.0

    children = Vector{BB.AbstractNode}()
    push!(children, child1)
    push!(children, child2)

    BB.push!(tree, children)

    @test tree.node_counter == 3
    @test length(tree.nodes) == 2
    @test length(tree.processed) == 1
end

next_node = BB.next_node(tree)

@testset "Child 1 with model" begin
    @test !isnothing(next_node.model)
    @test JuMP.lower_bound(next_node.reference_map[x]) == 0.
    @test JuMP.upper_bound(next_node.reference_map[x]) == 0.5
    @test next_node.dual_bound == 0.0
end

@testset "Bounding Child 1 node" begin
    BB.bound!(next_node)
    @test next_node.solution_status == MOI.OPTIMAL
    @test next_node.dual_bound == 0.0
    @test BB.node_solution(next_node, x) == 0.0
end
