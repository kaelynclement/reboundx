/**
 * @file    yorp_effect.c
 * @brief   Spin-up from YORP effect
 * @author  Kaelyn Clement <kaelynclement@gmail.com>
 *
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Radiation Forces$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 Kaelyn Clement
 * Implementation Paper    Clement et al., in prep.
 * Based on                `Veras and Scheeres., 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.2437V/abstract>`_.
 * Python Example          `YORPEffect.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/YORP_Effect.ipynb>`_.
 * ======================= ===============================================
 *
 * This calculates the change in rotation frequency (spin-up) of a body due to the YORP effect.
 * If the tensile strength is set, the user can also define a heartbeat function within problem.c
 * to determine how to handle particle fragmentation.
 * 
 * **Effect Parameters**
 *
 * ====================================== =========== ===========================================================
 * Field (C type)                         Required    Description
 * ====================================== =========== ===========================================================
 * yorp_solar_radiation_constant (float)  Yes         The solar radiation constant.
 * yorp_lstar (float)                     Yes         Luminosity of sim's star in units of solar luminosity.
 * ====================================== =========== ===========================================================
 *
 * **Particle Parameters**
 *
 * =============================== =========== ==================================================================
 * Field (C type)                  Required    Description
 * =============================== =========== ==================================================================
 * particles[i].r (float)          Yes         Physical radius of a body.
 * yorp_body_density (float)       Yes         Density of an object.
 * yorp_rotation_frequency (float) Yes         Rotation frequency of a spinning object.
 * o.a (float)                     Yes         Semimajor axis of body orbit.
 * o.e (float)                     Yes         Eccentricity of body orbit.
 * yorp_tensile_strength (float)   No          Tensile strength of a body.
 * yorp_obliquity (float)          No          Obliquity of an object in radians.
 * yorp_alpha (float)              No          Constant relating the z and epsilon YORP coefficients.
 * yorp_beta (float)               No          Fitting coefficient for obliquity-dependent YORP.
 * =============================== =========== ==================================================================
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "reboundx.h"

void rebx_yorp_effect(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){

    // only apply yorp effect to non-variational particles
    const int _N_real = sim->N - sim->N_var;

    // gravitational constant
    double G = sim->G;

	for(int i=1; i<_N_real; i++){
		struct reb_particle* p = &sim->particles[i];

        // orbital elements needed to calculate yorp effect
        struct reb_orbit o = reb_orbit_from_particle(G, *p, sim->particles[0]);
        const double a = o.a; // semi-major axis
        const double e = o.e; // eccentricity

        // particle radius needed to calculate yorp effect
        const double r = p->r;

        // yorp effect parameters needed for calculation
        const double* c_body = rebx_get_param(sim->extras, p->ap, "yorp_c_body");
        const double* phi = rebx_get_param(sim->extras, operator->ap, "yorp_solar_radiation_constant");
        const double* density = rebx_get_param(sim->extras, p->ap, "yorp_body_density");
        const double* lstar = rebx_get_param(sim->extras, operator->ap, "yorp_lstar");
        double* rotation_frequency = rebx_get_param(sim->extras, p->ap, "yorp_rotation_frequency");

        // parameters needed to include obliquity effects
        double* obliquity = rebx_get_param(sim->extras, p->ap, "yorp_obliquity");
        const double* alpha = rebx_get_param(sim->extras, p->ap, "yorp_alpha");
        const double* beta = rebx_get_param(sim->extras, p->ap, "yorp_beta");
        
        // check all parameters are defined before calculating yorp effect
        if (c_body != NULL && phi != NULL && density != NULL && lstar != NULL && rotation_frequency != NULL){

            const double frac = ((3.*(*c_body)*(*phi))/(4.*M_PI*(*density)*r*r*a*a*sqrt(1.-e*e)))*(*lstar);

            // check parameters needed to include the obliquity effects
            if (obliquity != NULL && alpha != NULL && beta != NULL){
                *rotation_frequency += frac*(cos(2.*(*obliquity))+*beta)*dt;
                *obliquity += (*alpha/(*rotation_frequency))*frac*sin(2.*(*obliquity))*dt;

                // keep obliquity between 0 and pi
                if (*obliquity > M_PI){
                    *obliquity -= M_PI;
                }
                if (*obliquity < 0.0){
                    *obliquity += M_PI;
                }
            }
            else {
                *rotation_frequency += frac*dt;
            }
        }
	}
}