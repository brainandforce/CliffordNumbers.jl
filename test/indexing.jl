@testset "BitIndex" begin
    import CliffordNumbers: signbit_of_square, nondegenerate_square, sign_of_square
    import CliffordNumbers: signbit_of_mult, nondegenerate_mult, sign_of_mult
    a = BitIndex(Val{VGA(3)}(), 1)
    b = BitIndex(Val{VGA(3)}(), 2)
    c = BitIndex(Val{VGA(3)}(), 3)
    @test BitIndex(Val{VGA(3)}(), 1, 2) === BitIndex{VGA(3)}(false, UInt(3))
    @test BitIndex(Val{VGA(3)}(), 2, 1) === BitIndex{VGA(3)}(true, UInt(3))
    @test BitIndex(Val{VGA(3)}(), 1, 2) === -BitIndex(Val{VGA(3)}(), 2, 1)
    @test abs(BitIndex(Val{VGA(3)}(), 1, 2)) === BitIndex(Val{VGA(3)}(), 1, 2)
    @test abs(BitIndex(Val{VGA(3)}(), 2, 1)) === BitIndex(Val{VGA(3)}(), 1, 2)
    # Sign manipulation
    @test copysign(+a, +1) === +a
    @test copysign(+a, -1) === -a
    @test copysign(-a, +1) === +a
    @test copysign(-a, -1) === -a
    @test flipsign(+b, +1) === +b
    @test flipsign(+b, -1) === -b
    @test flipsign(-b, +1) === -b
    @test flipsign(-b, -1) === +b
    @test copysign(+c, +a) === +c
    @test copysign(+c, -a) === -c
    @test copysign(-c, +a) === +c
    @test copysign(-c, -a) === -c
    @test flipsign(+c, +b) === +c
    @test flipsign(+c, -b) === -c
    @test flipsign(-c, +b) === -c
    @test flipsign(-c, -b) === +c
    # Euclidean multiplications
    @test signbit_of_mult(a, b) === false
    @test signbit_of_mult(b, a) === true
    @test signbit_of_mult(-a, b) === true
    @test signbit_of_mult(a, -b) === true
    @test signbit_of_mult(-a, -b) === false
    @test signbit_of_square(BitIndex(Val{VGA(3)}())) === false
    @test sign_of_square(BitIndex(Val{VGA(3)}())) > 0
    @test signbit_of_square(a) === false
    @test sign_of_square(a) > 0
    @test signbit_of_square(BitIndex(Val{VGA(3)}(), 1, 2)) === true
    @test sign_of_square(BitIndex(Val{VGA(3)}(), 1, 2)) < 0
    @test signbit_of_square(BitIndex(Val{VGA(3)}(), 1, 2, 3)) === true
    @test sign_of_square(BitIndex(Val{VGA(3)}(), 1, 2, 3)) < 0
    @test a * b === BitIndex(Val{VGA(3)}(), 1, 2)
    @test b * a === BitIndex(Val{VGA(3)}(), 2, 1)
    @test b * a === -BitIndex(Val{VGA(3)}(), 1, 2)
    @test CliffordNumbers.has_wedge(a, BitIndex(Val{VGA(3)}())) === true
    @test CliffordNumbers.has_wedge(a, b) === true
    @test CliffordNumbers.has_wedge(b, a) === true
    @test CliffordNumbers.has_wedge(a, a) === false
    @test CliffordNumbers.has_wedge(a, a*b) === false
    @test CliffordNumbers.has_wedge(a*b, a) === false
    @test CliffordNumbers.has_wedge(a*b, b) === false
    @test CliffordNumbers.has_wedge(a*b, a*b) === false
    @test CliffordNumbers.has_wedge(a*b, b*a) === false
    @test CliffordNumbers.has_wedge(a, b*c) === true
    @test CliffordNumbers.has_wedge(a*b, c) === true
    @test CliffordNumbers.has_wedge(a, b, c) === true
    @test CliffordNumbers.has_wedge(a, b, b) === false
    @test CliffordNumbers.has_wedge(a, b, c, b) === false
    @test nondegenerate_mult(a, b) === true
    @test nondegenerate_square(a*b) === true
    # Degenerate multiplications with Cl(3, 1, 2)
    QF = Signature(6, 0b001000, 0b110000, 1)
    @test nondegenerate_mult(BitIndex(Val{VGA(3)}(), 1, 3), BitIndex(Val{VGA(3)}(), 2, 3)) === true
    @test nondegenerate_mult(BitIndex(Val{QF}(), 5), BitIndex(Val{QF}(), 5)) === false
    @test nondegenerate_mult(BitIndex(Val{QF}(), 1, 5), BitIndex(Val{QF}(), 1, 5)) === false
    @test nondegenerate_mult(BitIndex(Val{QF}(), 1, 5), BitIndex(Val{QF}(), 2, 5)) === false
    @test nondegenerate_mult(BitIndex(Val{QF}(), 2, 5), BitIndex(Val{QF}(), 1, 5)) === false
    @test nondegenerate_mult(BitIndex(Val{QF}(), 6), BitIndex(Val{QF}(), 5)) === true
    @test nondegenerate_mult(BitIndex(Val{QF}(), 1, 6), BitIndex(Val{QF}(), 1, 5)) === true
    @test nondegenerate_mult(BitIndex(Val{QF}(), 1, 6), BitIndex(Val{QF}(), 2, 5)) === true
    @test nondegenerate_mult(BitIndex(Val{QF}(), 2, 6), BitIndex(Val{QF}(), 1, 5)) === true
    @test nondegenerate_square(BitIndex(Val{QF}(), 1, 3)) === true
    @test nondegenerate_square(BitIndex(Val{QF}(), 1, 5)) === false
    @test nondegenerate_square(BitIndex(Val{QF}(), 5, 6)) === false
    @test Base.Broadcast.broadcastable(BitIndex(Val(VGA(3)))) === tuple(BitIndex(Val(VGA(3))))
    @test eval(Meta.parse(repr(BitIndex(Val{VGA(3)}())))) === BitIndex(Val(VGA(3)))
    @test eval(Meta.parse(repr(BitIndex(Val{VGA(3)}(), 1, 2)))) === BitIndex(Val(VGA(3)), 1, 2)
    @test eval(Meta.parse(repr(BitIndex(Val{VGA(3)}(), 2, 1)))) === BitIndex(Val(VGA(3)), 2, 1)
end

@testset "BitIndices" begin
    # Tests to check that types don't proliferate like crazy
    @test BitIndices(CliffordNumber{VGA(3)}) === BitIndices{VGA(3),CliffordNumber{VGA(3)}}()
    @test BitIndices(CliffordNumber{VGA(3),Float32}) === BitIndices(CliffordNumber{VGA(3)})
    @test BitIndices(CliffordNumber{VGA(3),Complex{Int}}) === BitIndices(CliffordNumber{VGA(3)})
    @test BitIndices(CliffordNumber{VGA(3),Bool,8}) === BitIndices(CliffordNumber{VGA(3)})
    APS_bivector_indices = [
        BitIndex(Val{VGA(3)}(), 1, 2),
        BitIndex(Val{VGA(3)}(), 1, 3),
        BitIndex(Val{VGA(3)}(), 2, 3)
    ]
    @test BitIndices{VGA(3),KVector{2,VGA(3)}}() == APS_bivector_indices
    @test BitIndices(KVector{2,VGA(3)}(4,2,0)) == APS_bivector_indices
    @test BitIndices{VGA(3),KVector{2,VGA(3)}}() == BitIndices(KVector{2,VGA(3)}(4,2,0))
    @test grade.(BitIndices(VGA(3))) == count_ones.(0:7)
    @test scalar_index(zero(CliffordNumber{VGA(3)})) === BitIndex(Val{VGA(3)}())
    @test pseudoscalar_index(zero(CliffordNumber{VGA(3)})) === BitIndex(Val{VGA(3)}(), 1, 2, 3)
    @test all(map(-, BitIndices{VGA(3)}()) .== (-).(BitIndices{VGA(3)}()))
    @test Broadcast.BroadcastStyle(BitIndices) === Broadcast.Style{Tuple}()
    @test Broadcast.BroadcastStyle(TransformedBitIndices) === Broadcast.Style{Tuple}()
    # Efficient == methods
    @test BitIndices(KVector{1,VGA(3)}) == BitIndices{VGA(3), KVector{1,VGA(3)}}()
    @test BitIndices(KVector{1,VGA(3)}) == BitIndices{VGA(3), KVector{1,VGA(3),Int}}()
    @test BitIndices(KVector{1,VGA(3)}) == BitIndices{VGA(3), KVector{1,VGA(3),Int,3}}()
end

@testset "Transformed BitIndices" begin
    k = KVector{2,VGA(3)}(4, 2, 0)
    e = EvenCliffordNumber{VGA(3)}(k)
    x = CliffordNumber{VGA(3)}(k)
    @test all(reverse.(BitIndices(x)) .== Iterators.map(reverse, BitIndices(x)))
    @test all(reverse.(BitIndices(x)) .== ReversedBitIndices(x))
    @test all(reverse.(BitIndices(e)) .== Iterators.map(reverse, BitIndices(e)))
    @test all(reverse.(BitIndices(e)) .== ReversedBitIndices(e))
    @test all(reverse.(BitIndices(k)) .== Iterators.map(reverse, BitIndices(k)))
    @test all(reverse.(BitIndices(k)) .== ReversedBitIndices(k))
    @test all(grade_involution.(BitIndices(x)) .== Iterators.map(grade_involution, BitIndices(x)))
    @test all(grade_involution.(BitIndices(x)) .== GradeInvolutedBitIndices(x))
    @test all(grade_involution.(BitIndices(e)) .== Iterators.map(grade_involution, BitIndices(e)))
    @test all(grade_involution.(BitIndices(e)) .== GradeInvolutedBitIndices(e))
    @test all(grade_involution.(BitIndices(k)) .== Iterators.map(grade_involution, BitIndices(k)))
    @test all(grade_involution.(BitIndices(k)) .== GradeInvolutedBitIndices(k))
    @test all(conj.(BitIndices(x)) .== Iterators.map(conj, BitIndices(x)))
    @test all(conj.(BitIndices(x)) .== ConjugatedBitIndices(x))
    @test all(conj.(BitIndices(e)) .== Iterators.map(conj, BitIndices(e)))
    @test all(conj.(BitIndices(e)) .== ConjugatedBitIndices(e))
    @test all(conj.(BitIndices(k)) .== Iterators.map(conj, BitIndices(k)))
    @test all(conj.(BitIndices(k)) .== ConjugatedBitIndices(k))
end

@testset "Indexing" begin
    k = KVector{2,VGA(3)}(4, 2, 0)
    k_inf = KVector{1,STA}(1//0, 2, 3, 4)
    x = CliffordNumber{VGA(3)}(0, 0, 0, 4, 0, 2, 0, 0)
    @test iszero(k[BitIndex(Val{VGA(3)}())])
    @test iszero(k[BitIndex(Val{VGA(3)}(), 1)])
    @test iszero(k[BitIndex(Val{VGA(3)}(), 2)])
    @test iszero(k[BitIndex(Val{VGA(3)}(), 3)])
    @test iszero(k[BitIndex(Val{VGA(3)}(), 1, 2, 3)])
    @test k[BitIndices(k)] === k
    @test k[BitIndices(VGA(3))] === x
    @test x[BitIndices(k)] === k
    @test x[BitIndices(KVector{2,VGA(3)})] === k
    # Testing for cases where the first element is 1//0
    # This caused problems before commit 528636f4496785253be9807b17df1f028ef7a5f0
    @test k_inf[BitIndex(Val(STA))] === 0//1
    @test k_inf[BitIndex(Val(STA), 0)] === 1//0
    @test k_inf[BitIndex(Val(STA), 1)] === 2//1
end

@testset "Type lengths" begin
    @test nblades(Int) === 1
    @test nblades(ComplexF64) === 1
    @test nblades(KVector{2,STA}) === nblades(KVector{2,STA,Int,6})
    @test nblades(zero(KVector{2,STA})) === nblades(KVector{2,STA,Int,6})
    @test nblades(CliffordNumber{STA}) === nblades(CliffordNumber{STA,Int,16})
    @test nblades(zero(CliffordNumber{STA})) === nblades(CliffordNumber{STA,Int,16})
    @test nblades(EvenCliffordNumber{STA}) === nblades(EvenCliffordNumber{STA,Int,8})
    @test nblades(zero(EvenCliffordNumber{STA})) === nblades(EvenCliffordNumber{STA,Int,8})
end
