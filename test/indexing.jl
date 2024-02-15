@testset "BitIndex and BitIndices" begin
    @test BitIndex{VGA(3)}(1,2) === -BitIndex{VGA(3)}(2,1)
    @test abs(BitIndex{VGA(3)}(1,2)) === BitIndex{VGA(3)}(1,2)
    @test abs(BitIndex{VGA(3)}(2,1)) === BitIndex{VGA(3)}(1,2)
    aps_bivector_indices = [BitIndex{VGA(3)}(1,2), BitIndex{VGA(3)}(1,3), BitIndex{VGA(3)}(2,3)]
    @test BitIndices{VGA(3),KVector{2,VGA(3)}}() == aps_bivector_indices
    @test BitIndices(KVector{2,VGA(3)}(4,2,0)) == aps_bivector_indices
    @test BitIndices{VGA(3),KVector{2,VGA(3)}}() == BitIndices(KVector{2,VGA(3)}(4,2,0))
    @test grade.(BitIndices(VGA(3))) == count_ones.(0:7)
end

@testset "Indexing" begin
    k = KVector{2,VGA(3)}(4, 2, 0)
    x = CliffordNumber{VGA(3)}(0, 0, 0, 4, 0, 2, 0, 0)
    @test iszero(k[BitIndex{VGA(3)}()])
    @test iszero(k[BitIndex{VGA(3)}(1)])
    @test iszero(k[BitIndex{VGA(3)}(2)])
    @test iszero(k[BitIndex{VGA(3)}(3)])
    @test iszero(k[BitIndex{VGA(3)}(1,2,3)])
    @test k[BitIndices(k)] === k
    @test k[BitIndices(VGA(3))] === x
    @test x[BitIndices(k)] === k
    @test x[BitIndices(KVector{2,VGA(3)})] === k
end
