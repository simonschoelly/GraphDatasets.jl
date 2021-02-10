

@testset "TUDatasets" begin

    graphs = loadgraphs(TUDatasets.AIDSDataset())

    @test ng(graphs) == 2000
end
