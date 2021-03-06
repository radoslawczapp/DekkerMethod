
# DBL_EPSILON - jest to dokładność dla liczb zmiennoprzecinkowych
# dla Float64 jest to 2.2204460492503131e-16
DBL_EPSILON = 2.2204460492503131e-16

function DIFFSIGN(x::Float64, y::Float64)
    if (x <=0 && y >= 0) || (x >= 0 && y <= 0)
        return true
    else
        return false
    end
end

function fun(x::Float64)
    return 1.0 / (x - 3.0) - 6.0
end

function between(x::Float64, a::Float64, b::Float64)
    if b > a
        return (x >= a && x <= b)
    else
        return (x >= b && x <= a)
    end
end

function lfun(b::Float64, a::Float64, fb::Float64, fa::Float64)
    if fb != fa
        return b - fb*(b - a) / (fb - fa) # metoda siecznych
    elseif fa != 0
        return Inf
    else
        return b
    end
end

function hfun(b::Float64, c::Float64)
    if c > b
        return b + abs(b*DBL_EPSILON)
    else
        return b - abs(b*DBL_EPSILON)
    end
end

function mfun(b::Float64, c::Float64)
    return 0.5*(b + c) # metoda bisekcji
end

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

function DekkerA(x0::Float64, x1::Float64, Eps::Float64)
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
    fxp = fun(x0)
    fx = fun(x1)
    
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
    
    xk = x0
    fxk = fxp
    x = x1
    iter = 1
    
    while abs(b - c) > 2*Eps
        iter = iter + 1
        lambda = lfun(b, a, fb, fa)
        xp = x
        x = vfun(lambda, b, c)
        fxp = fx
        fx = fun(x)

        if DIFFSIGN(fxp, fx) == true
            xk = xp
            fxk = fxp
        end
        
        if abs(fx) <= abs(fxk)
            a = b
            b = x
            c = xk
            fa = fb
            fb = fx
        else
            b = xk
            a = c = x
            fa = fx
            fb = fxk
        end       
    end
    println("Number of iterations: ", iter)
    return b
end

println("f(x) = 1/(x-3)-6\n")
result = DekkerA(3.01, 4.0, 1e-12)
println("x0 = ", result)
