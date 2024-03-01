@testset "Hamming weight tools" begin
    import CliffordNumbers:
        isevil, isodious, evil_number, odious_number, next_of_hamming_weight, hamming_number
    # Sequences pulled from OEIS
    first_16_evil = [0, 3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30]
    first_16_odious = [1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28, 31]
    hamming_weight_3 = [7, 11, 13, 14, 19, 21, 22, 25, 26, 28, 35, 37, 38, 41, 42, 44]
    @test all(isevil, first_16_evil)
    @test all(isodious, first_16_odious)
    @test evil_number.(1:16) == first_16_evil
    @test odious_number.(1:16) == first_16_odious
    @test all(x -> next_of_hamming_weight(hamming_weight_3[x]) == hamming_weight_3[x+1], 1:15)
    @test all(x -> count_ones(x) == 3, hamming_weight_3)
    @test hamming_number.(3, 1:16) == hamming_weight_3
    @test hamming_number.(1, 1:8) == [2^(n-1) for n in 1:8]
    @test all(x -> isone(count_ones(x)), hamming_number.(1, 1:8))
end

@testset "Type relations" begin
    @test CliffordNumber <: AbstractCliffordNumber
    @test EvenCliffordNumber <: AbstractCliffordNumber
    @test OddCliffordNumber <: AbstractCliffordNumber
    @test KVector <: AbstractCliffordNumber
end

@testset "Represented grades" begin
    @test nonzero_grades(KVector{2,APS}) === 2:2
    @test nonzero_grades(zero(KVector{2,APS})) === 2:2
    @test nonzero_grades(CliffordNumber{APS}) === 0:3
    @test nonzero_grades(zero(CliffordNumber{APS})) === 0:3
    @test RepresentedGrades(KVector{2,APS})[0:3] == [false, false, true, false]
    @test RepresentedGrades(CliffordNumber{APS})[0:3] == trues(4)
end

@testset "Printing/display" begin
    import CliffordNumbers.short_typename
    @test short_typename(zero(CliffordNumber{APS,Float64,8})) === CliffordNumber{APS,Float64}
    @test short_typename(zero(CliffordNumber{APS,Float64})) === CliffordNumber{APS,Float64}
    @test short_typename(zero(EvenCliffordNumber{APS,Int,4})) === EvenCliffordNumber{APS,Int}
    @test short_typename(zero(OddCliffordNumber{APS,Int})) === OddCliffordNumber{APS,Int}
    @test short_typename(zero(KVector{1,APS,Bool,3})) === KVector{1,APS,Bool}
    @test short_typename(zero(KVector{2,APS,Bool})) === KVector{2,APS,Bool}
end
