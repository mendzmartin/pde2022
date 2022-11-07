# algunas funciones auxiliares

@inline function mysign(a)
    return (1.0.*(a .> 0.0) + (-1.0).* (a .< 0.0))
end

@inline minmod(x, y) = 0.5*(mysign(x)+mysign(y))*min(abs(x),abs(y))
@inline function minmod(a, b, c)
    sgnbc = (mysign(b)+mysign(c)) #Esto es 2 si ambos son positivos, -2 si ambos son negativos, 0 en otro caso
    sgnac = (mysign(a)+mysign(c)) #Esto es 2 si ambos son positivos, -2 si ambos son negativos, 0 en otro caso
    return 0.25*sgnbc*abs(sgnac)*min(abs(a), abs(b),abs(c))
end

@inline function minmod(w, x, y, z) 
    0.125*(mysign(w) + mysign(x))*abs((mysign(w)+mysign(y))*(mysign(w)+mysign(z)))*min(abs(w),abs(x),abs(y),abs(z))
end


#método de Kurganov-Tadmor

function KT!(dfields, fields, par, t)
    #Los parámetros son h, θ, funciones auxiliares y vectores auxiliares
    eqpars, h::Float64, θ::Float64, Fx!, MaxSpeed, N::Int64, N_FIELDS::Int32, auxvectors = par
    Dm, D, Dp, u_mm, u_mp, u_pm, u_pp, F_mm, F_mp, F_pm, F_pp, H_m, H_p = auxvectors

    for idx in 1:N
        idxll = mod(((idx-2) - 1),N) + 1
        idxl = mod(((idx-1) - 1),N) + 1
        idxr = mod(((idx+1) - 1),N) + 1
        idxrr = mod(((idx+2) - 1),N) + 1
        
        fll = @view fields[idxll,:]
        fl = @view fields[idxl,:]
        f = @view fields[idx,:]
        fr = @view fields[idxr,:]
        frr = @view fields[idxrr,:]
        
        @. Dm = minmod(0.5 *(f - fll), θ*(f-fl), θ*(fl-fll))    # *(h)
        @. D  = minmod(0.5 *(fr - fl), θ*(fr-f), θ*(f-fl))      # *(h)
        @. Dp = minmod(0.5*(frr-f), θ*(frr-fr), θ*(fr-f))       # *(h)
        @. u_mm = fl + 0.5*Dm   # *(1/h)
        @. u_mp = f - 0.5*D     # *(1/h)
        @. u_pm = f + 0.5*D     # *(1/h)
        @. u_pp = fr - 0.5*Dp   # *(1/h)
        
        Fx!(F_mm, u_mm, eqpars)
        Fx!(F_mp, u_mp, eqpars)
        Fx!(F_pm, u_pm, eqpars)
        Fx!(F_pp, u_pp, eqpars)
        
        
        a_m::Float64 = max(MaxSpeed(u_mm, eqpars), MaxSpeed(u_mp, eqpars))       
        a_p::Float64 = max(MaxSpeed(u_pm, eqpars), MaxSpeed(u_pp, eqpars))        

        
        @. H_m = 0.5 * (F_mp + F_mm) - 0.5 * a_m * (u_mp - u_mm)
        @. H_p = 0.5 * (F_pp + F_pm) - 0.5 * a_p * (u_pp - u_pm)
        
        @. dfields[idx, :] = -h*(H_p - H_m)
    end
end

function createKTauxvectors(N_FIELDS)
    D = Array{Float64}(undef, N_FIELDS)
    Dm = copy(D)
    Dp = copy(D)
    umm = copy(D)
    ump = copy(D)
    upm = copy(D)
    upp = copy(D)
    Fmm = copy(D)
    Fmp = copy(D)
    Fpm = copy(D)
    Fpp = copy(D)
    Hm = copy(D)
    Hp = copy(D)
    return (Dm, D, Dp, umm, ump, upm, upp, Fmm, Fmp, Fpm, Fpp, Hm, Hp)
end

#==================================MP5=====================================#

#Reconstrucción
function MP5reconstruction!(Vl, Vjmm, Vjm, Vj, Vjp, Vjpp, N_Fields)
    B1 = 0.0166666666666666667  # = 1/60
    B2 = 1.3333333333333333333  # = 4/3
    eps = 1e-10
    ALPHA = 4.0
    #=Vjmm = V[1]
    Vjm = V[2]
    Vj = V[3]
    Vjp = V[4]
    Vjpp = V[5]=#
    for i in 1:N_Fields
        # This is the original interpolation. All that follows is the application of limiters to treat shocks
        Vor = B1*(2.0*Vjmm[i] - 13.0*Vjm[i] + 47.0*Vj[i] + 27*Vjp[i] - 3.0*Vjpp[i])
        # mp = monotonicity preserving. It's the median between v_j, v_(j+1) and an upper limit v^UL = v_j+ALPHA(v_j-v_(j-1))
        Vmp = Vj[i] + minmod(Vjp[i]-Vj[i], ALPHA*(Vj[i]-Vjm[i]))
        if ((Vor-Vj[i])*(Vor-Vmp)) < eps    # this condition is equivalent to asking vl in [vj, v^{MP}]
            Vl[i] = Vor                     # vl = v^{L}_{j+1/2}
        else
            djm1 = Vjmm[i] - 2.0*Vjm[i] + Vj[i]
            dj = Vjm[i] - 2*Vj[i] + Vjp[i]
            djp1 = Vj[i] - 2.0*Vjp[i] + Vjpp[i]
            dm4jph = minmod(4*dj - djp1, 4*djp1-dj, dj, djp1)  # ph = plus half (+1/2)
            dm4jmh = minmod(4*dj - djm1, 4*djm1-dj, dj, djm1)  # mh = minus half (-1/2)
            # d^{M4}_{j+1/2} = \minmod(4d_{j}-d_{j+1},4d_{j+1}-d_{j}, d_{j}, d_{j+1})
            Vul = Vj[i] + ALPHA*(Vj[i] - Vjm[i])    # upper limit
            Vav = 0.5*(Vj[i] + Vjp[i])              # average
            Vmd = Vav - 0.5*dm4jph                  # Vmedian
            Vlc = Vj[i] + 0.5*(Vj[i]-Vjm[i]) + B2*dm4jmh
            Vmin = max(min(Vj[i], Vjp[i], Vmd), min(Vj[i], Vul, Vlc));
            Vmax = min(max(Vj[i], Vjp[i], Vmd), max(Vj[i], Vul, Vlc));
            Vl[i] = Vor + minmod(Vmin-Vor, Vmax-Vor) #this places Vor between Vmin and Vmax
        end
    end
end 

#Implementación de MP5 con el Flux Splitting de Lax

function mp5!(dfields, fields, par, t) # j is the grid position
    #asumimos u unidimensional por ahora
    par_eq, h, N, N_Fields, Fx!, Speed_max, auxvectors = par
    F_Mm3, F_Mm2, F_Mm1, F_M, F_Mp1, F_Mp2, F_Mp3, F_Pm3, F_Pm2, F_Pm1, F_P, F_Pp1, F_Pp2, F_Pp3, F_LP, F_LM, F_RP, F_RM, H_m, H_p = auxvectors
    

    #nota: f minuscula o u se usa para hablar de campos, F mayúscula para hablar de Flujos.
    
    for idx in 1:N
        #first we defined shifted indices
        idxm3 = mod(((idx-3) - 1),N) + 1
        idxm2 = mod(((idx-2) - 1),N) + 1
        idxm1 = mod(((idx-1) - 1),N) + 1
        idxp1 = mod(((idx+1) - 1),N) + 1
        idxp2 = mod(((idx+2) - 1),N) + 1
        idxp3 = mod(((idx+3) - 1),N) + 1
        
    
        um3 = @view fields[idxm3,:]
        um2 = @view fields[idxm2,:]
        um1 = @view fields[idxm1,:]
        u   = @view fields[idx,:]
        up1 = @view fields[idxp1,:]
        up2 = @view fields[idxp2,:]
        up3 = @view fields[idxp3,:]
        
        S_MAX = max(Speed_max(up3, par_eq), Speed_max(um3, par_eq), 
            Speed_max(up2, par_eq), Speed_max(um2, par_eq), Speed_max(up1, par_eq), 
            Speed_max(um1, par_eq), Speed_max(u, par_eq)) #maximum speed
        
        Fx!(F_Pm3, um3, par_eq)
        Fx!(F_Pm2, um2, par_eq)
        Fx!(F_Pm1, um1, par_eq)
        Fx!(F_P, u, par_eq)
        Fx!(F_Pp1, up1, par_eq)
        Fx!(F_Pp2, up2, par_eq)
        Fx!(F_Pp3, up3, par_eq)
        
        
        @. F_Mm3 = 0.5 * (F_Pm3 - S_MAX * um3)
        @. F_Mm2 = 0.5 * (F_Pm2 - S_MAX * um2)
        @. F_Mm1 = 0.5 * (F_Pm1 - S_MAX * um1)
        @. F_M   = 0.5 * (F_P   - S_MAX * u)
        @. F_Mp1 = 0.5 * (F_Pp1 - S_MAX * up1)
        @. F_Mp2 = 0.5 * (F_Pp2 - S_MAX * up2)
        @. F_Mp3 = 0.5 * (F_Pp3 - S_MAX * up3)
        @. F_Pm3 = 0.5 * (F_Pm3 + S_MAX * um3)
        @. F_Pm2 = 0.5 * (F_Pm2 + S_MAX * um2)
        @. F_Pm1 = 0.5 * (F_Pm1 + S_MAX * um1)
        @. F_P   = 0.5 * (F_P   + S_MAX * u)
        @. F_Pp1 = 0.5 * (F_Pp1 + S_MAX * up1)
        @. F_Pp2 = 0.5 * (F_Pp2 + S_MAX * up2)
        @. F_Pp3 = 0.5 * (F_Pp3 + S_MAX * up3)
    
        MP5reconstruction!(F_RM, F_Mp2, F_Mp1,  F_M,  F_Mm1, F_Mm2, N_Fields)
        MP5reconstruction!(F_LM, F_Pm3, F_Pm2, F_Pm1, F_P,  F_Pp1, N_Fields)
        MP5reconstruction!(F_LP, F_Pm2, F_Pm1,  F_P,  F_Pp1, F_Pp2, N_Fields)
        MP5reconstruction!(F_RP, F_Mp3, F_Mp2, F_Mp1, F_M,  F_Mm1, N_Fields)
        
        @. H_p = F_LP + F_RP
        @. H_m = F_LM + F_RM
        
        @. dfields[idx, :] = -h*(H_p - H_m)
        
    end
    
end

function createMP5auxvectors(N_FIELDS)
    F_P = Array{Float64}(undef, N_FIELDS)
    F_P = copy(F_P)
    F_M = copy(F_P)
    F_Pm3 = copy(F_P)
    F_Pm2 = copy(F_P)
    F_Pm1 = copy(F_P)
    F_Pp1 = copy(F_P)
    F_Pp2 = copy(F_P)
    F_Pp3 = copy(F_P)
    F_Mm3 = copy(F_P)
    F_Mm2 = copy(F_P)
    F_Mm1 = copy(F_P)
    F_Mp1 = copy(F_P)
    F_Mp2 = copy(F_P)
    F_Mp3 = copy(F_P)

    F_LP = copy(F_P)
    F_LM = copy(F_P)
    F_RM = copy(F_P)
    F_RP = copy(F_P)
    H_m = copy(F_P)
    H_p = copy(F_P)
    return (F_Mm3, F_Mm2, F_Mm1, F_M, F_Mp1, F_Mp2, F_Mp3, F_Pm3, F_Pm2, F_Pm1, F_P, F_Pp1, F_Pp2, F_Pp3, F_LP, F_LM, F_RP, F_RM, H_m, H_p)
end