@testset "BitIndex" begin
    a = BitIndex(VGA(3), 1)
    b = BitIndex(VGA(3), 2)
    c = BitIndex(VGA(3), 3)
    @test BitIndex(VGA(3), 1, 2) === BitIndex{VGA(3)}(false, UInt(3))
    @test BitIndex(VGA(3), 2, 1) === BitIndex{VGA(3)}(true, UInt(3))
    @test BitIndex(VGA(3), 1, 2) === -BitIndex(VGA(3), 2, 1)
    @test abs(BitIndex(VGA(3), 1, 2)) === BitIndex(VGA(3), 1, 2)
    @test abs(BitIndex(VGA(3), 2, 1)) === BitIndex(VGA(3), 1, 2)
    # Euclidean multiplications
    @test CliffordNumbers.signbit_of_mult(a, b) === false
    @test CliffordNumbers.signbit_of_mult(b, a) === true
    @test CliffordNumbers.signbit_of_mult(-a, b) === true
    @test CliffordNumbers.signbit_of_mult(a, -b) === true
    @test CliffordNumbers.signbit_of_mult(-a, -b) === false
    @test CliffordNumbers.signbit_of_square(BitIndex(VGA(3))) === false
    @test CliffordNumbers.sign_of_square(BitIndex(VGA(3))) > 0
    @test CliffordNumbers.signbit_of_square(a) === false
    @test CliffordNumbers.sign_of_square(a) > 0
    @test CliffordNumbers.signbit_of_square(BitIndex(VGA(3), 1, 2)) === true
    @test CliffordNumbers.sign_of_square(BitIndex(VGA(3), 1, 2)) < 0
    @test CliffordNumbers.signbit_of_square(BitIndex(VGA(3), 1, 2, 3)) === true
    @test CliffordNumbers.sign_of_square(BitIndex(VGA(3), 1, 2, 3)) < 0
    @test a * b === BitIndex(VGA(3), 1, 2)
    @test b * a === BitIndex(VGA(3), 2, 1)
    @test b * a === -BitIndex(VGA(3), 1, 2)
    @test CliffordNumbers.has_wedge(a, BitIndex(VGA(3))) === true
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
    @test CliffordNumbers.nondegenerate_mult(a, b) === true
    @test CliffordNumbers.nondegenerate_square(a*b) === true
    # Degenerate multiplications
    QF = QuadraticForm{3,1,2}
    @test CliffordNumbers.nondegenerate_mult(BitIndex(APS, 1, 3), BitIndex(APS, 2, 3)) === true
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 5), BitIndex(QF, 5)) === false
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 1, 5), BitIndex(QF, 1, 5)) === false
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 1, 5), BitIndex(QF, 2, 5)) === false
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 2, 5), BitIndex(QF, 1, 5)) === false
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 6), BitIndex(QF, 5)) === true
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 1, 6), BitIndex(QF, 1, 5)) === true
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 1, 6), BitIndex(QF, 2, 5)) === true
    @test CliffordNumbers.nondegenerate_mult(BitIndex(QF, 2, 6), BitIndex(QF, 1, 5)) === true
    @test CliffordNumbers.nondegenerate_square(BitIndex(QF, 1, 3)) === true
    @test CliffordNumbers.nondegenerate_square(BitIndex(QF, 1, 5)) === false
    @test CliffordNumbers.nondegenerate_square(BitIndex(QF, 5, 6)) === false
end

@testset "BitIndices" begin
    aps_bivector_indices = [BitIndex(VGA(3), 1, 2), BitIndex(VGA(3), 1, 3), BitIndex(VGA(3), 2, 3)]
    @test BitIndices{VGA(3),KVector{2,VGA(3)}}() == aps_bivector_indices
    @test BitIndices(KVector{2,VGA(3)}(4,2,0)) == aps_bivector_indices
    @test BitIndices{VGA(3),KVector{2,VGA(3)}}() == BitIndices(KVector{2,VGA(3)}(4,2,0))
    @test grade.(BitIndices(VGA(3))) == count_ones.(0:7)
    @test scalar_index(zero(CliffordNumber{VGA(3)})) === BitIndex(VGA(3))
    @test pseudoscalar_index(zero(CliffordNumber{VGA(3)})) === BitIndex(VGA(3), 1, 2, 3)
    @test all(map(-, BitIndices{VGA(3)}()) .== (-).(BitIndices{VGA(3)}()))
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
    x = CliffordNumber{VGA(3)}(0, 0, 0, 4, 0, 2, 0, 0)
    @test iszero(k[BitIndex(VGA(3))])
    @test iszero(k[BitIndex(VGA(3), 1)])
    @test iszero(k[BitIndex(VGA(3), 2)])
    @test iszero(k[BitIndex(VGA(3), 3)])
    @test iszero(k[BitIndex(VGA(3), 1, 2, 3)])
    @test k[BitIndices(k)] === k
    @test k[BitIndices(VGA(3))] === x
    @test x[BitIndices(k)] === k
    @test x[BitIndices(KVector{2,VGA(3)})] === k
end
