
# DBL_EPSILON - jest to dokładność dla liczb zmiennoprzecinkowych
# dla Float64 jest to 2.2204460492503131e-16
DBL_EPSILON = 2.2204460492503131e-16

function DIFFSIGN(x, y)
    if (x <=0 && y >= 0) || (x >= 0 && y <= 0)
        return true
    else
        return false
    end
end

function fun(x)
#     return 1 / (x - 3) - 6
    return (x + 3)*((x - 1)^2)
end

function between(x, a, b)
    if b > a
        return (x >= a && x <= b)
    else
        return (x >= b && x <= a)
    end
end

function lfun(b, a, fb, fa)
    if fb != fa
        return b - fb*(b - a) / (fb - fa) # metoda siecznych
    elseif fa != 0
        return Inf
    else
        return b
    end
end

function hfun(b, c)
    if c > b
        return b + abs(b*DBL_EPSILON)
    else
        return b - abs(b*DBL_EPSILON)
    end
end

function mfun(b, c)
    return 0.5*(b + c) # metoda bisekcji
end
"""
function vfun(l::Float64, b::Float64, c::Float64)
    h = hfun(b, c)
    m = mfun(b, c)
    
    if between(l, h, m) == true
        return l
    elseif abs(l - b) <= abs(b*DBL_EPSILON)
        return h
    else
        return m
    end
end
"""

function wfun(l, b, c)
    h = hfun(b, c)
    m = mfun(b, c)
    
    if between(l, h, m) == true
         return l
    elseif (fun(abs(l - b)) <= abs(b*DBL_EPSILON)) && (between(l, b, m) == false)
        return h
    else
        return m
    end
end

function ffun(a, b, fa, fb)
     return (fa - fb) / (a - b)
end

function rfun(b, a, d, fb, fa, fd)
    alpha = ffun(b, d, fb, fd)*fa
    beta = ffun(a, d, fa, fd)*fb
    
    if beta != alpha
        return b - beta*(b - a) / (beta - alpha)
    elseif alpha != 0
        return Inf
    else
        return 0 # beta == alpha == 0
    end
end

function DekkerM(x0, x1, Eps)
    """
    b - ostatnie przyblizenie wyniku
    c - kontrapunkt b, punkt w którym funkcja f ma przeciwny znak niż w punkcie b
    a – poprzednia wartość punktu a, używana do wyliczania następnego punktu metodą siecznych
    
    Metoda bisekcji z punktów b i c tworzy punkt m pomiędzy nimi na środku przedziału.
    
    Wyliczany jest ciąg xi, którego ostatni element oznaczany jest przez x, a poprzedni przez xp.
    
    xk - ostatni punkt w ciągu, który ma różny znak niż x.
    
    Punkt x wyliczany jest dwoma metodami:
        - siecznych
        - bisekcji
    i wybierany jest ten obliczony z metody siecznych jeśli leży pomiędzy punktem b
    (ze względów dokładnościowych z pewną poprawką) a punktem m wyliczonym z bisekcji.

    Jeżeli f(x) czy f(xk) leży bliżej zera i jeśli f(x) leży bliżej zera
    wtedy b ma wartość x, c ma wartość xk, w przeciwnym razie zamiana.
    """
    d = NaN
    fd = NaN
    
    fxp = fun(x0)
    fx = fun(x1)
    
    if x0 == x1
        return fx
    end
    
    if DIFFSIGN(fx, fxp) == false
        return 0
    end
    
    if abs(fx) <= abs(fxp)
        b = x1
        a = c = x0
        fa = fxp
        fb = fx
    else
        b = x0
        a = c = x1
        fa = fx
        fb = fxp
    end
    
    xk = xp = x0
    fxk = fxp
    x = x1
    iter = 1
    age = 0
    bp = b
    cp = c
    ap = a
    
    while abs(b - c) > 2*Eps
        iter = iter + 1
        age = age + 1
        
        if abs(b - c) <= (0.5 + 2 * DBL_EPSILON)*(abs(bp - cp) + abs(b*DBL_EPSILON))
            age = 1
        end
        
        xp = x
        
        if age <= 2
            lambda = lfun(b, a, fb, fa)
            
            if abs(lambda - b) < abs(b*DBL_EPSILON)
                break
            end
            x = wfun(lambda, b, c)
        elseif age == 3
            rho = rfun(b, a, d, fb, fa, fd)
            
            if abs(rho - b) < abs(b*DBL_EPSILON)
                break
            end
            x = wfun(rho, b, c)
        else
            x = mfun(b, c)
        end
        
        fxp = fx
        fx = fun(x)

        if DIFFSIGN(fxp, fx) == true
            xk = xp
            fxk = fxp
        end
        
        bp = b
        fbp = fb
        ap = a
        fap = fa
        cp = c
        
        if abs(fx) <= abs(fxk)
            a = b
            fa = fb
            b = x
            fb = fx
            c = xk
        else
            b = xk
            fb = fxk
            a = c = x
            fa = fx
        end
        
        if (b == x) || (b == bp)
            d = ap
            fd = fap
        else
            d = ap
            fd = fdp
        end
    end
    println("Number of iterations: ", iter)
    return b
end

println("f(X) = (x+3)(x-1)^2")
result = DekkerM(-4, 4/3, 1e-12)
# println("f(x) = 1/(x-3)-6\n")
# result = DekkerM(3.01, 4, 1e-12)
println("x0 = ", result)
