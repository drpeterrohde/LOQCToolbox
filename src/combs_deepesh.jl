using Combinatorics, Symbolics, SymbolicUtils, LinearAlgebra, SymPy

@vars t

function SumPermInt(termsLeft, termsRight)
    perms = collect(permutations(termsRight))
    terms = length(termsLeft)

    expanded_arr = []
    for i in perms
        push!(expanded_arr, termsLeft*i)
    end
    #println(expanded_arr)

    sum_of_arr = sum(expanded_arr)
    #println(sum_of_arr[1])

    final_term = integrate(sum_of_arr[1], (t,-Inf,Inf))
    #println(final_term)

    return final_term
end

function SumPermInt(terms)
    sum = SumPermInt(conj(terms), terms)
    return(sum)
end

# Example
v=[exp(-t^2) exp(-t^2)]
dummy = SumPermInt(v,v)
println(dummy)
