@testset "Random number generation" begin
    n_samples = 8192
    rng1 = Random.Xoshiro(420)
    randn_VGA2 = [randn(rng1, KVector{1,VGA(2)}) for _ in 1:n_samples]
    randn_VGA3 = [randn(rng1, KVector{1,VGA(3)}) for _ in 1:n_samples]
    moduli_VGA2 = abs2.(randn_VGA2)
    moduli_VGA3 = abs2.(randn_VGA3)
    moduli_VGA2_mean = mean(moduli_VGA2)
    moduli_VGA3_mean = mean(moduli_VGA3)
    moduli_VGA2_var = var(moduli_VGA2, mean=moduli_VGA2_mean)
    moduli_VGA3_var = var(moduli_VGA3, mean=moduli_VGA3_mean)
    # Check that the moduli match a chi-squared distribution using a Z-test
    # Mean should be equal to the number of dimensions
    # Variance should be twice the mean
    @test abs(moduli_VGA2_mean - 2) / sqrt(moduli_VGA3_var + 2*2) < 2
    @test abs(moduli_VGA3_mean - 3) / sqrt(moduli_VGA3_var + 2*3) < 2
end
