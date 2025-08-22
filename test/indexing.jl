@testset "BladeIndex" begin
    import CliffordNumbers: signbit_of_square, nondegenerate_square, sign_of_square
    import CliffordNumbers: signbit_of_mult, nondegenerate_mult, sign_of_mult
    a = BladeIndex(Val{VGA(3)}(), 1)
    b = BladeIndex(Val{VGA(3)}(), 2)
    c = BladeIndex(Val{VGA(3)}(), 3)
    @test BladeIndex(Val{VGA(3)}(), 1, 2) === BladeIndex{VGA(3)}(false, UInt(3))
    @test BladeIndex(Val{VGA(3)}(), 2, 1) === BladeIndex{VGA(3)}(true, UInt(3))
    @test BladeIndex(Val{VGA(3)}(), 1, 2) === -BladeIndex(Val{VGA(3)}(), 2, 1)
    @test BladeIndex(Val{VGA(2)}(), 1, 2) != BladeIndex(Val{VGA(3)}(), 1, 2)
    @test abs(BladeIndex(Val{VGA(3)}(), 1, 2)) === BladeIndex(Val{VGA(3)}(), 1, 2)
    @test abs(BladeIndex(Val{VGA(3)}(), 2, 1)) === BladeIndex(Val{VGA(3)}(), 1, 2)
    @test signature(a) === VGA(3)
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
    # With other numeric types
    @test flipsign(1, a) === 1
    @test flipsign(1, -a) === -1
    @test flipsign(-1, a) === -1
    @test flipsign(-1, -a) === 1
    @test copysign(2.0, b) === 2.0
    @test copysign(-2.0, b) === 2.0
    @test copysign(2.0, -b) === -2.0
    @test copysign(-2.0, b) === -2.0
    # Euclidean multiplications
    @test signbit_of_mult(a, b) === false
    @test signbit_of_mult(b, a) === true
    @test signbit_of_mult(-a, b) === true
    @test signbit_of_mult(a, -b) === true
    @test signbit_of_mult(-a, -b) === false
    @test signbit_of_square(BladeIndex(Val{VGA(3)}())) === false
    @test sign_of_square(BladeIndex(Val{VGA(3)}())) > 0
    @test signbit_of_square(a) === false
    @test sign_of_square(a) > 0
    @test signbit_of_square(BladeIndex(Val{VGA(3)}(), 1, 2)) === true
    @test sign_of_square(BladeIndex(Val{VGA(3)}(), 1, 2)) < 0
    @test signbit_of_square(BladeIndex(Val{VGA(3)}(), 1, 2, 3)) === true
    @test sign_of_square(BladeIndex(Val{VGA(3)}(), 1, 2, 3)) < 0
    @test a * b === BladeIndex(Val{VGA(3)}(), 1, 2)
    @test b * a === BladeIndex(Val{VGA(3)}(), 2, 1)
    @test b * a === -BladeIndex(Val{VGA(3)}(), 1, 2)
    @test adjoint(a) === reverse(a)
    @test (a * b)' === b * a
    @test CliffordNumbers.has_wedge(a, BladeIndex(Val{VGA(3)}())) === true
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
    @test nondegenerate_mult(BladeIndex(Val{VGA(3)}(), 1, 3), BladeIndex(Val{VGA(3)}(), 2, 3))
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 5), BladeIndex(Val{QF}(), 5)) === false
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 1, 5), BladeIndex(Val{QF}(), 1, 5)) === false
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 1, 5), BladeIndex(Val{QF}(), 2, 5)) === false
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 2, 5), BladeIndex(Val{QF}(), 1, 5)) === false
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 6), BladeIndex(Val{QF}(), 5)) === true
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 1, 6), BladeIndex(Val{QF}(), 1, 5)) === true
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 1, 6), BladeIndex(Val{QF}(), 2, 5)) === true
    @test nondegenerate_mult(BladeIndex(Val{QF}(), 2, 6), BladeIndex(Val{QF}(), 1, 5)) === true
    @test nondegenerate_square(BladeIndex(Val{QF}(), 1, 3)) === true
    @test nondegenerate_square(BladeIndex(Val{QF}(), 1, 5)) === false
    @test nondegenerate_square(BladeIndex(Val{QF}(), 5, 6)) === false
    @test Base.Broadcast.broadcastable(BladeIndex(Val(VGA(3)))) === tuple(BladeIndex(Val(VGA(3))))
end

@testset "BladeIndices" begin
    # Tests to check that types don't proliferate like crazy
    @test BladeIndices(CliffordNumber{VGA(3)}) === BladeIndices{VGA(3),CliffordNumber{VGA(3)}}()
    @test BladeIndices(CliffordNumber{VGA(3),Float32}) === BladeIndices(CliffordNumber{VGA(3)})
    @test BladeIndices(CliffordNumber{VGA(3),Complex{Int}}) === BladeIndices(CliffordNumber{VGA(3)})
    @test BladeIndices(CliffordNumber{VGA(3),Bool,8}) === BladeIndices(CliffordNumber{VGA(3)})
    APS_bivector_indices = [
        BladeIndex(Val{VGA(3)}(), 1, 2),
        BladeIndex(Val{VGA(3)}(), 1, 3),
        BladeIndex(Val{VGA(3)}(), 2, 3)
    ]
    @test BladeIndices{VGA(3),KVector{2,VGA(3)}}() == APS_bivector_indices
    @test BladeIndices{VGA(3)}(KVector{2,VGA(3)}(4,2,0)) == APS_bivector_indices
    @test BladeIndices(KVector{2,VGA(3)}(4,2,0)) == APS_bivector_indices
    @test BladeIndices{VGA(3),KVector{2,VGA(3)}}() == BladeIndices(KVector{2,VGA(3)}(4,2,0))
    @test grade.(BladeIndices(VGA(3))) === count_ones.(NTuple{8,Int}(0:7))
    @test scalar_index(zero(CliffordNumber{VGA(3)})) === BladeIndex(Val{VGA(3)}())
    @test pseudoscalar_index(zero(CliffordNumber{VGA(3)})) === BladeIndex(Val{VGA(3)}(), 1, 2, 3)
    @test all(map(-, BladeIndices{VGA(3)}()) .== (-).(BladeIndices{VGA(3)}()))
    @test Broadcast.BroadcastStyle(BladeIndices) === Broadcast.Style{Tuple}()
    @test Broadcast.BroadcastStyle(CliffordNumbers.CGNBladeIndices) === Broadcast.Style{Tuple}()
    # Efficient == methods
    @test BladeIndices(KVector{1,VGA(3)}) == BladeIndices{VGA(3), KVector{1,VGA(3)}}()
    @test BladeIndices(KVector{1,VGA(3)}) == BladeIndices{VGA(3), KVector{1,VGA(3),Int}}()
    @test BladeIndices(KVector{1,VGA(3)}) == BladeIndices{VGA(3), KVector{1,VGA(3),Int,3}}()
end

@testset "CGNBladeIndices" begin
    import CliffordNumbers: CGNBladeIndices, CyclicGradeNegation
    k = KVector{2,VGA(3)}(4, 2, 0)
    e = EvenCliffordNumber(k)
    x = CliffordNumber(k)
    Tk = CliffordNumbers.blade_indices_type(k)
    Te = CliffordNumbers.blade_indices_type(e)
    Tx = CliffordNumbers.blade_indices_type(x)
    @test all(identity.(BladeIndices(x)) .== Iterators.map(identity, BladeIndices(x)))
    @test identity.(BladeIndices(x)) isa BladeIndices{VGA(3),Tx}
    @test all(identity.(BladeIndices(e)) .== Iterators.map(identity, BladeIndices(e)))
    @test identity.(BladeIndices(e)) isa BladeIndices{VGA(3),Te}
    @test all(identity.(BladeIndices(k)) .== Iterators.map(identity, BladeIndices(k)))
    @test identity.(BladeIndices(k)) isa BladeIndices{VGA(3),Tk}
    @test all(reverse.(BladeIndices(x)) .== Iterators.map(reverse, BladeIndices(x)))
    @test reverse.(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,false,false,true}
    @test all(reverse.(BladeIndices(e)) .== Iterators.map(reverse, BladeIndices(e)))
    @test reverse.(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,false,false,true}
    @test all(reverse.(BladeIndices(k)) .== Iterators.map(reverse, BladeIndices(k)))
    @test reverse.(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,false,false,true}
    @test all(grade_involution.(BladeIndices(x)) .== Iterators.map(grade_involution, BladeIndices(x)))
    @test grade_involution.(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,false,true,false}
    @test all(grade_involution.(BladeIndices(e)) .== Iterators.map(grade_involution, BladeIndices(e)))
    @test grade_involution.(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,false,true,false}
    @test all(grade_involution.(BladeIndices(k)) .== Iterators.map(grade_involution, BladeIndices(k)))
    @test grade_involution.(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,false,true,false}
    @test all(conj.(BladeIndices(x)) .== Iterators.map(conj, BladeIndices(x)))
    @test conj.(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,false,true,true}
    @test all(conj.(BladeIndices(e)) .== Iterators.map(conj, BladeIndices(e)))
    @test conj.(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,false,true,true}
    @test all(conj.(BladeIndices(k)) .== Iterators.map(conj, BladeIndices(k)))
    @test conj.(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,false,true,true}
    @test all((-).(BladeIndices(x)) .== Iterators.map(-, BladeIndices(x)))
    @test (-).(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,true,false,false}
    @test all((-).(BladeIndices(e)) .== Iterators.map(-, BladeIndices(e)))
    @test (-).(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,true,false,false}
    @test all((-).(BladeIndices(k)) .== Iterators.map(-, BladeIndices(k)))
    @test (-).(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,true,false,false}
    @test all(-reverse.(BladeIndices(x)) .== Iterators.map(reverse, -BladeIndices(x)))
    @test -reverse.(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,true,false,true}
    @test all(-reverse.(BladeIndices(e)) .== Iterators.map(reverse, -BladeIndices(e)))
    @test -reverse.(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,true,false,true}
    @test all(-reverse.(BladeIndices(k)) .== Iterators.map(reverse, -BladeIndices(k)))
    @test -reverse.(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,true,false,true}
    @test all(-grade_involution.(BladeIndices(x)) .== Iterators.map(grade_involution, -BladeIndices(x)))
    @test -grade_involution.(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,true,true,false}
    @test all(-grade_involution.(BladeIndices(e)) .== Iterators.map(grade_involution, -BladeIndices(e)))
    @test -grade_involution.(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,true,true,false}
    @test all(-grade_involution.(BladeIndices(k)) .== Iterators.map(grade_involution, -BladeIndices(k)))
    @test -grade_involution.(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,true,true,false}
    @test all(-conj.(BladeIndices(x)) .== Iterators.map(conj, -BladeIndices(x)))
    @test -conj.(BladeIndices(x)) isa CGNBladeIndices{VGA(3),Tx,true,true,true}
    @test all(-conj.(BladeIndices(e)) .== Iterators.map(conj, -BladeIndices(e)))
    @test -conj.(BladeIndices(e)) isa CGNBladeIndices{VGA(3),Te,true,true,true}
    @test all(-conj.(BladeIndices(k)) .== Iterators.map(conj, -BladeIndices(k)))
    @test -conj.(BladeIndices(k)) isa CGNBladeIndices{VGA(3),Tk,true,true,true}
    # Map definitions that don't generate Tuple
    @test map(+, BladeIndices(k)) === CGNBladeIndices{VGA(3),Tk,false,false,false}()
    @test map(identity, BladeIndices(k)) === CGNBladeIndices{VGA(3),Tk,false,false,false}()
    @test map(-, BladeIndices(k)) === CGNBladeIndices{VGA(3),Tk,true,false,false}()
    @test map(grade_involution, BladeIndices(k)) === CGNBladeIndices{VGA(3),Tk,false,true,false}()
    @test map(reverse, BladeIndices(k)) === CGNBladeIndices{VGA(3),Tk,false,false,true}()
    @test map(conj, BladeIndices(k)) === CGNBladeIndices{VGA(3),Tk,false,true,true}()
    @test +BladeIndices(e) === BladeIndices(e)
    # Not sure if we'll use this type on its own, honestly
    @test CyclicGradeNegation{true,false,false}()(k) === -k
    @test CyclicGradeNegation{false,true,false}()(k) === k
    @test CyclicGradeNegation{false,false,true}()(k) === -k
end

@testset "Indexing" begin
    k = KVector{2,VGA(3)}(4, 2, 0)
    k_inf = KVector{1,STA}(1//0, 2, 3, 4)
    x = CliffordNumber{VGA(3)}(0, 0, 0, 4, 0, 2, 0, 0)
    @test iszero(k[BladeIndex(Val{VGA(3)}())])
    @test iszero(k[BladeIndex(Val{VGA(3)}(), 1)])
    @test iszero(k[BladeIndex(Val{VGA(3)}(), 2)])
    @test iszero(k[BladeIndex(Val{VGA(3)}(), 3)])
    @test iszero(k[BladeIndex(Val{VGA(3)}(), 1, 2, 3)])
    @test k[BladeIndices(k)] === k
    @test k[BladeIndices(VGA(3))] === x
    @test x[BladeIndices(k)] === k
    @test x[BladeIndices(KVector{2,VGA(3)})] === k
    # Testing for cases where the first element is 1//0
    # This caused problems before commit 528636f4496785253be9807b17df1f028ef7a5f0
    @test k_inf[BladeIndex(Val(STA))] === 0//1
    @test k_inf[BladeIndex(Val(STA), 0)] === 1//0
    @test k_inf[BladeIndex(Val(STA), 1)] === 2//1
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

@testset "Custom printing" begin
    @test eval(Meta.parse(repr(BladeIndex(Val{VGA(3)}())))) === BladeIndex(Val(VGA(3)))
    @test eval(Meta.parse(repr(BladeIndex(Val{VGA(3)}(), 1, 2)))) === BladeIndex(Val(VGA(3)), 1, 2)
    @test eval(Meta.parse(repr(BladeIndex(Val{VGA(3)}(), 2, 1)))) === BladeIndex(Val(VGA(3)), 2, 1)
    inds = BladeIndices(CliffordNumber{STA})
    @test eval(Meta.parse(repr(inds))) === inds
    @test eval(Meta.parse(repr(grade_involution.(inds)))) === grade_involution.(inds)
    @test eval(Meta.parse(repr(reverse.(inds)))) === reverse.(inds)
    @test eval(Meta.parse(repr(conj.(inds)))) === conj.(inds)
    @test eval(Meta.parse(repr(-inds))) === -inds
    @test eval(Meta.parse(repr(-grade_involution.(inds)))) === -grade_involution.(inds)
    @test eval(Meta.parse(repr(-adjoint.(inds)))) === -adjoint.(inds)
    @test eval(Meta.parse(repr(-conj.(inds)))) === -conj.(inds)
end
