module Transformations

export ri_dq, dq_ri

"""
    ri_dq(δ::T) where {T<:Number}

RI to DQ transformation matrix.
Uses the reference frame of Kundur page 852.
"""
function ri_dq(δ::T) where {T<:Number}
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

"""
    dq_ri(δ::T) where {T<:Number}

DQ to RI transformation matrix.
Uses the reference frame of Kundur page 852.
"""
function dq_ri(δ::T) where {T<:Number}
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

end # module