@testset "ChainRulesCore extension" begin
    k = KVector{1,VGA(3)}(4, 2, 0)
    l = KVector{2,VGA(3)}(0, 6, 9)
    test_rrule(*, k, l)
end
