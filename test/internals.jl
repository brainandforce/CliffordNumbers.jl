@testset "Internals" begin
    @test CliffordNumbers.intlog2(3//4) === round(Int, log2(3//4))
    @test CliffordNumbers.intlog2(0.75) === round(Int, log2(0.75))
    @test CliffordNumbers.intlog2(1) === 0
    @test CliffordNumbers.intlog2(1//1) === 0
    @test CliffordNumbers.intlog2(1.0) === 0
    @test CliffordNumbers.intlog2(3//2) === round(Int, log2(3//2))
    @test CliffordNumbers.intlog2(1.5) === round(Int, log2(1.5))
    # Should always pass unless we specialize on integers
    @test all(CliffordNumbers.intlog2(n) === CliffordNumbers.intlog2(Float64(n)) for n in 0:65535)
end

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
