using MTH229
using Plots

function bisec(f::Function, x1::Number, x2::Number, e::AbstractFloat=1e-7, max::Integer=100)
    fx1 = f(x1)
    fx2 = f(x2)
    fx1*fx2 <= 0 || error("The functions do not have a real root solution.")
    iter = 0
    local valor
    while x2 - x1 > e
        iter += 1
        iter != max 
        valor = (x2+x1)/2
        fvalor = f(valor)
        if fvalor == 0 
            break
        elseif fx1*fvalor > 0 
            x1 = valor
            fx1 = fvalor
        else 
            x2 = valor
        end
    end
    return valor
end

function falsi(f::Function, x1::Number, x2::Number, 
               maxiter::Integer=50, ey=2eps(Float64), ex::AbstractFloat=1e-7)
    fx1 = f(x1)
    fx2 = f(x2)
    fx1*fx2 <= 0 || error("The functions do not have a real root solution.")
    iter = 0
    while iter < maxiter
        iter += 1
        xk = x2-(fx2*((x2-x1)/(fx2-fx1)))
        fxk = f(xk)
        if min(abs(xk-x1), abs(xk-x2)) < ex
            return xk
        end
        if abs(fxk) < ey
            return xk
        end 
        if sign(fx1*fxk) == 1
            x1 = xk
            fx1 = fxk
        else 
            x2 = xk
            fx2 = fxk
        end
    end
    @sprintf("xk = %f",xk)
    @sprintf("fxk = %f",fxk)
end  