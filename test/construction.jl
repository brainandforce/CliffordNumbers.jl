@testset "Zero elements" begin
    @test zero(KVector{0,APS}) == KVector{0,APS,Bool}(0)
    @test zero(KVector{1,APS}) == KVector{1,APS,Bool}(0, 0, 0)
    @test zero(KVector{2,APS}) == KVector{2,APS,Bool}(0, 0, 0)
    @test zero(KVector{3,APS}) == KVector{3,APS,Bool}(0)
    @test zero(EvenCliffordNumber{APS}) == EvenCliffordNumber{APS,Bool}(0, 0, 0, 0)
    @test zero(OddCliffordNumber{APS}) == OddCliffordNumber{APS,Bool}(0, 0, 0, 0)
    @test zero(CliffordNumber{APS}) == CliffordNumber{APS,Bool}(0, 0, 0, 0, 0, 0, 0, 0)
end

@testset "Units" begin
    # one
    @test one(KVector{0,APS}) == KVector{0,APS,Bool}(1)
    @test one(KVector{1,APS}) == KVector{0,APS,Bool}(1)
    @test one(KVector{2,APS}) == KVector{0,APS,Bool}(1)
    @test one(KVector{3,APS}) == KVector{0,APS,Bool}(1)
    @test one(EvenCliffordNumber{APS}) == EvenCliffordNumber{APS,Bool}(1, 0, 0, 0)
    @test one(OddCliffordNumber{APS}) == KVector{0,APS,Bool}(1)
    @test one(CliffordNumber{APS}) == CliffordNumber{APS,Bool}(1, 0, 0, 0, 0, 0, 0, 0)
    # oneunit
    @test oneunit(KVector{0,APS}) == KVector{0,APS,Bool}(1)
    @test oneunit(EvenCliffordNumber{APS}) == EvenCliffordNumber{APS,Bool}(1, 0, 0, 0)
    @test oneunit(CliffordNumber{APS}) == CliffordNumber{APS,Bool}(1, 0, 0, 0, 0, 0, 0, 0)
    # these should throw since they can't represent 1
    @test_throws InexactError oneunit(KVector{1,APS})
    @test_throws InexactError oneunit(KVector{2,APS})
    @test_throws InexactError oneunit(KVector{3,APS})
    @test_throws InexactError oneunit(OddCliffordNumber{APS})
end
