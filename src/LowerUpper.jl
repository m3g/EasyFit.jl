@with_kw struct lower
    a = nothing
    b = nothing
    c = nothing
    d = nothing
end

function lower(a::Number, b::Number, c::Number)
    lower(a,b,c,nothing)
end

function lower(a::Number, b::Number)
    lower(a,b,nothing,nothing)
end

function lower(a::Number)
    lower(a,nothing,nothing,nothing)
end

@with_kw struct upper
    a = nothing
    b = nothing
    c = nothing
    d = nothing
end

function upper(a::Number, b::Number, c::Number)
    upper(a,b,c,nothing)
end

function upper(a::Number, b::Number)
    upper(a,b,nothing,nothing)
end

function upper(a::Number)
    upper(a,nothing,nothing,nothing)
end

export lower, upper
