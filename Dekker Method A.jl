
DBL_EPSILON = 2.2204460492503131e-16

function DIFFSIGN(x, y)
    if (x <=0 && y >= 0) || (x >= 0 && y <= 0)
        return true
    else
        return false
    end
end

function fun(x)
    return (1 / (x - 3)) - 6
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
        return b - fb*(b - a) / (fb - fa) # secant method
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
    return 0.5*(b + c) # bisection method
end

function vfun(l, b, c)
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

function DekkerA(x0, x1, Eps)    
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
    
    while abs(b - c) > 2*Eps
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
    return b
end

result = DekkerA(3.01, 4, 1e-12)
println("x0 = ", result)
