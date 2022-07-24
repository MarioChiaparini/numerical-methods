using MTH229
using Plots

#################
##ROOTS FINDING##
#################

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

function secant(f::Function, x1::Number, x2::Number, maxiter::Integer=100)
    fx1 = f(x1)
    fx2 = f(x2)
    
    if fx2*fx1 >= 0
        return nothing
    end 
    iter = 0
    
    while iter < maxiter
        
        iter += 1
        xk = x1 - f(x1)*(x2 - x1)/(f(x1) - f(x2))
        fxk = f(xk)
        
        if fx1*fxk < 0
            x1 = x1
            x2 = xk  
        
        elseif fx2*fxk < 0
            x1 = xk
            x2 = x2
        
        elseif fxk == 0
            return xk 
        
        else 
            return nothing
        end 
    end
    return x1 - f(x1)*(x2 - x1)/(f(x2) - f(x1))
end

function newrap(f::Function, x0::Number, maxiter::Integer, tol::AbstractFloat=1e-5) 
    iter = 0
    df = diff(f)
    while abs(f(x0)) > tol || iter > maxiter
        iter += 1 
        xk =  x0-(f(x0)/df(x0))
        @sprintf("xk = %f",xk)
        @sprintf("fxk = %f",f(xk))
        x0 = xk
        return xk
    end
end

###########################
###NUMERICAL INTEGRATION###
###########################

#################
##EDO SOLUTION###
#################
#function rk4(f::Function, xn::AbstractFloat, x0::Number, y0::Number, n::Integer)
#    h = (xn-x0)/n
#    @sprintf("solution")
#    i=0
#    while i < n
#        i += 1
#        k1 = h * (f(x0,y0))
#        k2 = h * (f((x0+h/2), (y0+k1/2)))
#        k3 = h * (f((x0+h/2), (y0+k2/2)))
#        k4 = h * (f((x0+h), (y0+k3)))
#        k =  (k1+2*k2+2*k3+k4)/6
#        yn = y0 + k
#        @sprintf("--------------------------")
#        @sprintf("%.3f\t%.3f\t%.3f", x0,y0,yn)
#        y0 = yn
#        x0 = x0+h
#    @sprintf("\n x=%.3f, y=%3.f",xn,yn)
#    end
#    return xn,yn 
#end

#function rungekutta4(f, y0, t)
#    n = length(t)
#    y = zeros((n, length(y0)))
#    y[1,:] = y0
#    for i in 1:n-1
#        h = t[i+1] - t[i]
#        k1 = f(y[i,:], t[i])
#        k2 = f(k1 * h/2 .+ y[i,:], h/2 .+ y[i,:])
#        k3 = f(k2 * h/2 .+ y[i,:], h/2 .+ y[i,:])
#        k4 = f(k3 * h .+ y[i,:], h .+ y[i,:])
#        y[i+1,:] = (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4) .+ y[i,:]
#    end
#    return y
#end

function rgk4(dxdy, tn0, tnf, y0, n)
    tn = LinRange(tn0, tnf, n+1)
    y = zeros(n, length(y0))
    y[1,:] = y0
    for i in 1:n-1
        j = i + 1
        h = tn[j] - tn[i]
        k1 = f(y[i,:], t[i])
        k2 = f(k1 * h/2 .+ y[i,:], h/2 .+ y[i,:])
        k3 = f(k2 * h/2 .+ y[i,:], h/2 .+ y[i,:])
        k4 = f(k3 * h .+ y[i,:], h .+ y[i,:])
        y[j,:] = (h/6) * (k1 .+ 2*k2 .+ 2*k3 .+ k4) .+ y[i,:]
    end
    return y
end 