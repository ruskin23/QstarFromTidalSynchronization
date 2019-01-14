BEGIN {
    printf("age,R_star,L_star,M_star,R_tachocline,M_ini,M_conv,M_rad,I_conv,I_rad\n");
    sigma = 5.670373e-08 #Stefan-Boltzman constant in SI
    L_sun = 3.846e+26 #Solar luminosity in Watts
    R_sun = 695508000.0 #Solar radius in meters
    M_star = 0.06 #Stellar mass in solar masses.
    pi = atan2(0, -1)
}

{
    if($1 != "#") {
        age = 10.0**$1
        T_eff = $3;
        R_star = $2;
        L_star = (4.0 * pi * R_star**2 * R_sun**2 sigma * T_eff**4\
                  /\
                  L_sun);
        printf(\
            "%.1f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f\n",\
            age,\
            R_star,\
            L_star,\
            M_star,\
            0.0,\
            M_star,\
            M_star,\
            0.0,\
            M_star * R_star**2,\
            0.0\
        );
    }
}

END {
        age = 1e11
        printf(\
            "%.1f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f\n",\
            age,\
            R_star,\
            L_star,\
            M_star,\
            0.01,\
            M_star,\
            M_star,\
            0.01,\
            M_star * R_star**2,\
            0.01\
        );
}
