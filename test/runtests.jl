using CSV, DataFrames
using LinearAlgebra, Distributions
using Random
using SGCRNAs
using Test

num_genes=500; num_samples=20; num_modules=4;
Random.seed!(42)
genes_per_module = div(num_genes, num_modules)
gene_names = ["Gene_$(i)" for i in 1:num_genes]

ρ = 0.8
Σ = [ρ^abs(i-j) for i in 1:num_samples, j in 1:num_samples]
dist = MvNormal(zeros(num_samples), Symmetric(Σ))

expression_data = Array{Float64}(undef, num_genes, num_samples)
for module_id in 0:num_modules-1
    base_signal = rand(dist)
    for i in 1:genes_per_module
        gene_index = module_id * genes_per_module + i
        noise = rand(Normal(0, 0.1), num_samples)
        expression_data[gene_index, :] = base_signal + noise
    end
end

@testset "Function tests" begin
    CorData, GradData = (0, 0);
    @test_nowarn begin
        CorData, GradData = CGM(gene_names, expression_data);
    end
    
    clust, pos, edge_data = (0, 0, 0);
    redirect_stderr(devnull) do
        clust, pos, edge_data = SpectralClustering(CorData, GradData);
    end
    @test_nowarn isa(clust, Vector)
    
    d = 1; k = maximum(clust[d]);
    nw, new_pos, cnctdf, new_clust, score = (0, 0, 0, 0, 0);
    @test_nowarn begin
        nw, new_pos, cnctdf, new_clust, score = SetNetwork(edge_data, clust[d], pos, il=collect(1:k));
    end
    
    @test_nowarn begin
        DrawNetwork("./AllNetWork-0.5.png", nw, new_pos, cnctdf, new_clust, k, node_scores=score, edge_mode=:ALL, edge_threshold=0.5)
    end
    
    @test_nowarn begin
        CorPhenMod(DataFrame(expression_data, :auto), DataFrame(Matrix(expression_data'), :auto)[:,[3,1,4]], new_clust, "./CorPhenMod.png")
    end
end
