#
# Function to set bounds with optional user input 
#
function setbounds(vars, l, u, oney)
    # Throw error if the user set a bound to a variable which is not in list for this function
    for field in fieldnames(typeof(l))
        if !isnothing(getfield(l, field)) || !isnothing(getfield(u, field))
            found = false
            for var in vars
                if field == var.field
                    found = true
                    if var.type == Nothing
                        throw(ArgumentError("Bounds to intercept $(var.field) are not supported. For a constant use $(var.field) = const."))
                    end
                end
            end
            if !found
                throw(ArgumentError("A bound was set to variable '$field', but '$field' is not a variable of the current fit function."))
            end
        end
    end
    # Total number variables
    n = 0
    for var in vars
        n += var.dim
    end
    lower = Vector{typeof(ustrip(oney))}(undef, n)
    upper = Vector{typeof(ustrip(oney))}(undef, n)
    idim = 1
    for var in vars
        ltmp = check_lower_bound_input(var, getfield(l, var.field), oney)
        utmp = check_upper_bound_input(var, getfield(u, var.field), oney)
        @. lower[idim:idim+var.dim-1] = ltmp
        @. upper[idim:idim+var.dim-1] = utmp
        idim += var.dim
    end
    i = 0
    for var in vars
        i += 1
        if lower[i] > upper[i]
            throw(ArgumentError(" Error in bounds. Lower bound of '$(var.field)' greater than upper bound. "))
        end
    end
    return lower, upper
end

function check_lower_bound_input(var, value, oney)
    T = typeof(ustrip(oney))
    if isnothing(value)
        return [typemin(T) for _ in 1:var.dim]
    end
    if unit(eltype(value)) != unit(1)
        throw(ArgumentError("Currently the lower and upper bounds must be provided without units."))
    end
    if !(value isa var.type)
        throw(ArgumentError("Lower bound of $(var.field) must be of type $(var.type), got $(typeof(value))"))
    end
    if (value isa Vector) && (length(value) != var.dim)
        throw(ArgumentError("Lower bound of $(var.field) must be of dimension $(var.dim), got $(length(value))"))
    end
    return value
end

function check_upper_bound_input(var, value, oney)
    T = typeof(ustrip(oney))
    if isnothing(value)
        return [typemax(T) for _ in 1:var.dim]
    end
    if unit(eltype(value)) != unit(1)
        throw(ArgumentError("Currently the lower and upper bounds must be provided without units."))
    end
    if !(value isa var.type)
        throw(ArgumentError("Upper bound of $(var.field) must be of type $(var.type), got $(typeof(value))"))
    end
    if (value isa Vector) && (length(value) != var.dim)
        trhow(ArgumentError("Upper bound of $(var.field) must be of dimension $(var.dim), got $(length(value))"))
    end
    return value
end

@testitem "check_bounds" begin
    using EasyFit
    using Unitful
    x = rand(10); y = rand(10)
    @test_throws ArgumentError fitlinear(x,y; l=lower(a=3.0), u=upper(a=1.))
    @test_throws ArgumentError fitlinear(x,y; l=lower(a=3.0u"m"))
    @test_throws ArgumentError fitlinear(x,y; l=lower(b=3.0u"m"))
    @test_throws ArgumentError fitlinear(x,y; l=lower(c=3.0u"m"))
    @test_throws ArgumentError fitlinear(x,y; l=lower(a=[1.0, 2.0]))
    @test_throws ArgumentError fitlinear(x,y; u=upper(a=3.0u"m"))
    @test_throws ArgumentError fitlinear(x,y; u=upper(b=3.0u"m"))
    @test_throws ArgumentError fitlinear(x,y; u=upper(c=3.0u"m"))
    @test_throws ArgumentError fitlinear(x,y; u=upper(a=[1.0, 2.0]))
end




