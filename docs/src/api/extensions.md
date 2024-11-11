# Extensions

## [Quaternions.jl]

```@docs
(::Type{H})(::AbstractCliffordNumber{VGA(3)}) where H<:Quaternion
Quaternions.slerp(::AbstractCliffordNumber{VGA(3)}, ::AbstractCliffordNumber{VGA(3)}, ::Real)
Quaternions.slerp(::AbstractCliffordNumber{VGA(3)}, ::Quaternion, ::Real)
```

[Quaternions.jl]:   https://github.com/JuliaGeometry/Quaternions.jl
