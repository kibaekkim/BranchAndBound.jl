using Test
using BranchAndBound
using JuMP, GLPK
using MathOptInterface

const BB = BranchAndBound
const MOI = MathOptInterface

model = Model(GLPK.Optimizer)
@variable(model, 0 <= x <= 1)
@objective(model, Min, x)

tree = BB.AbstractTree()
node = BB.JuMPNode{BB.VariableBranch}(model)

@testset "Initial Tree" begin
    @test tree.node_counter == 0
    @test length(tree.nodes) == 0
    @test length(tree.processed) == 0
end

@testset "Initial Node" begin
    @test node.id == -1
    @test isnothing(node.parent)
    @test node.depth == 0
    @test isnothing(node.branch)
    @test node.solution_status == MOI.OPTIMIZE_NOT_CALLED
    @test isinf(-node.bound)
    @test !isnothing(node.model)
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
    @test JuMP.lower_bound(x) == 0.
    @test JuMP.upper_bound(x) == 1.
end

@testset "Bounding root node" begin
    BB.bound!(current_node)
    @test current_node.solution_status == MOI.OPTIMAL
    @test current_node.bound == 0.0
    @test JuMP.value(x) == 0.0
end

@testset "Marked root node as processed" begin
    BB.processed!(tree, current_node)
    @test isempty(tree.nodes)
    @test length(tree.processed) == 1
    @test tree.processed[1] == current_node
end

@testset "Branching at root" begin
    child1 = BB.create_child_node_with_lb(current_node, x, 0.5)
    # child1 = BB.JuMPNode(current_node)
    # @test child1.parent == current_node
    # @test child1.depth == 1
    # @test child1.bound == current_node.bound
    # child1.branch = BB.VariableBranch(Dict(x=>0.5),Dict{JuMP.VariableRef,Real}())

    child2 = BB.create_child_node_with_ub(current_node, x, 0.5)
    # child2 = BB.JuMPNode(current_node)
    # @test child2.parent == current_node
    # @test child2.depth == 1
    # @test child2.bound == current_node.bound
    # child2.branch = BB.VariableBranch(Dict{JuMP.VariableRef,Real}(),Dict(x=>0.5))

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
    print(next_node.model)
    @test JuMP.lower_bound(x) == 0.
    @test JuMP.upper_bound(x) == 0.5
    @test next_node.bound == 0.0
end

@testset "Bounding Child 1 node" begin
    BB.bound!(next_node)
    @test next_node.solution_status == MOI.OPTIMAL
    @test next_node.bound == 0.0
    @test JuMP.value(x) == 0.0
end
