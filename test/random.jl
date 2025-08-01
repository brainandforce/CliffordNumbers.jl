@testset "Random number generation" begin
    rng1 = Random.Xoshiro(420)
    randn_VGA2 = [randn(rng1, KVector{1,VGA(2)}) for _ in 1:1024]
    randn_VGA3 = [randn(rng1, KVector{1,VGA(3)}) for _ in 1:1024]
    moduli_VGA2 = abs2.(randn_VGA2)
    moduli_VGA3 = abs2.(randn_VGA3)
    moduli_VGA2_mean = mean(moduli_VGA2)
    moduli_VGA3_mean = mean(moduli_VGA3)
    moduli_VGA2_var = var(moduli_VGA2, mean=moduli_VGA2_mean)
    moduli_VGA3_var = var(moduli_VGA3, mean=moduli_VGA3_mean)
    # Check that the squares are part of the chi-squared distribution
end
