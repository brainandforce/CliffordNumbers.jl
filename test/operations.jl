@testset "Equality" begin
    x = CliffordNumber{APS,Float64}(1, 2, 3, 4, 5, 6, 7, 8)
    @test x == convert(CliffordNumber{APS,Int}, x)
    @test x !== convert(CliffordNumber{APS,Int}, x)
end

@testset "Addition" begin
    x = CliffordNumber{APS}(1, 2, 3, 4, 5, 6, 7, 8)
    y = CliffordNumber{APS,Float64}(9, -10, 11, -12, 13, 14, -15, 16)
    # Type promotion test
    @test x + y isa CliffordNumber{APS,Float64}
    # Equality test (mixed types)
    @test x + y == CliffordNumber{APS}(10, -8, 14, -8, 18, 20, -8, 24)
end

@testset "Geometric product" begin
    x = CliffordNumber{APS}(0, 2, 0, 0, 0, 0, 0, 0)
    y = CliffordNumber{APS}(0, 3, 4, 0, 0, 0, 0, 0)
    @test x * y == CliffordNumber{APS}(6, 0, 0, 8, 0, 0, 0, 0)
    # When multiplying vectors, the reverse should have a negative bivector
    @test y * x == CliffordNumber{APS}(6, 0, 0, -8, 0, 0, 0, 0)
end

@testset "Wedge product" begin
    x = CliffordNumber{APS}(0, 2, 0, 0, 0, 0, 0, 0)
    y = CliffordNumber{APS}(0, 3, 4, 0, 0, 0, 0, 0)
    five = 5 * one(CliffordNumber{APS})
    # Self wedges should be zero
    @test iszero(x ∧ x)
    @test iszero(y ∧ y)
    # Reversing the order should flip the sign for vectors
    @test x ∧ y == -(y ∧ x)
    # Ensure the behavior of scalars (both CliffordNumber and normal) are correct
    @test 5 ∧ x == 5 * x
    @test five ∧ x == x ∧ five
    @test 5 ∧ 5 == 25
end
