using Combinatorics, Symbolics, SymbolicUtils, QuadGK

@variables t

function SumPermInt(termsLeft, termsRight)
    sum = 0.0
    perms = collect(permutations(termsRight))
    terms = length(termsLeft)
    
    for perm in perms
        prod = 1.0
        for i = 1:terms
            f = eval(build_function(termsLeft[i]*perm[i],t))
          #  println(f(1.0))
            #println(quadgk(thisf,-10.0,10.0)[1])
            prod *= quadgk(f,-10,10)[1] # TODO: make range infinite
        end
        sum += prod
    end

    return(sum)
end

function SumPermInt(terms)
    sum = SumPermInt(conj(terms), terms)
    return(sum)
end

# TEST
v=[exp(-t^2) exp(-t^2)]