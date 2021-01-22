export convertF_C, convertF_K, convertin_m, convertm_in, convertPa_psi, convertpsi_Pa


"""
convertin_m(x) intakes a distance in inches and outputs that distance in meters.
"""
function convertin_m(x)
    return x*0.0254
end

"""
convertm_in(x) intakes a distance in meters and outputs that distance in inches.
"""
function convertm_in(x)
    return x/0.0254
end

function convertF_C(F)
    #Convert from Fahrenheit to Celcius
    return 5*(F-32)/9
end

function convertF_K(F)
    #Convert from Fahrenheit to Kelvin
    return convertF_C(F)+273.15
end

function convertpsi_Pa(psi)
    return psi*6894.7572931783
end

function convertPa_psi(Pa)
    return Pa/6894.7572931783
end