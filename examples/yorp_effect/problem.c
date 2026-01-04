/**
 * YORP Effect on a Single Asteroid
 * 
 * This example shows the effect of the YORP radiation force on the rotation frequency and obliquity of a single
 * asteroid in orbit around a red giant. The fragmentation process is defined in a custom heartbeat function and produces
 * a binary asteroid with a few smaller fragments called "shards". Since we are using test particles (where the mass
 * is set to 0), this is a simplified model and does not include binary interactions between the asteroid fragments. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);
int tsteps; // number of timesteps that have passed

int main(int argc, char* argv[]) {

    tsteps = 0; // initialize number of timesteps

    struct reb_simulation* sim = reb_simulation_create(); // create simulation

    sim->G = 4*M_PI*M_PI; // use units of AU, yr and solar masses
    sim->dt = 1.; // timestep for simulation in yrs
    sim->integrator = REB_INTEGRATOR_WHFAST; // integrator for sim
    sim->heartbeat = heartbeat; // function pointer for heartbeat

    // add star with mass of Sun to sim
    struct reb_particle star = {0};
    star.m = 1.;
    reb_simulation_add(sim, star);

    // asteroid orbital parameters
    double m = 0;
    double a = 5.; // semi-major axis of 5 AU
    double e = 0; // no eccentricity
    double inc = 0;
    double Omega = 0;
    double omega = 0;
    double f = 0;

    // create asteroid and add to sim
    struct reb_particle asteroid_1 = reb_particle_from_orbit(sim->G, star, m, a, e, inc, Omega, omega, f);
    reb_simulation_add(sim,asteroid_1);
    reb_simulation_move_to_com(sim);
    
    // pointers needed to set parameters and add YORP effect to sim
    struct reb_particle* const particles = sim->particles; // pointer for the particles in the sim
    struct rebx_extras* rebx = rebx_attach(sim); // attach REBOUNDx to the sim
    struct rebx_operator* yorp = rebx_load_operator(rebx, "yorp_effect"); // pointer for yorp operator
    
    // constant conversions
    double au_conv = 1.495978707e11; // meters in one AU
    double msun_conv = 1.9885e30; // kg in one solar mass
    double yr_conv = 31557600.0; // seconds in one year

    // effect parameters
    double phi = 1.e17/au_conv/msun_conv*yr_conv*yr_conv; // solar radiation constant in units of (solar mass AU yr^-2)
    double lstar = 1000.; // luminosity MUST be in units of solar luminosity

    // asteroid parameters
    double rotation_frequency = 0.0*yr_conv; // initial rotation frequency in rad/yr
    double c_body = 1./10./10.; // YORP coefficient (asymmetry) of the body
    double density = (2000.0*au_conv*au_conv*au_conv)/msun_conv; // density in units of solar mass AU^-3
    double radius = 100./au_conv; // radius in units of AU
    double sigma = 1.e3/msun_conv*au_conv*yr_conv*yr_conv; // tensile strength in units of solar mass AU^-1 yr^-2
    
    // asteroid obliquity parameters
    double obliquity = M_PI/6.; // initial obliquity in radians
    double alpha = 2.0/3.0; // constant relating the axial and obliquity YORP coefficients
    double beta = 1.0/3.0; // fitting coefficient

    // set parameters for the YORP effect
    rebx_set_param_double(rebx, &yorp->ap, "yorp_solar_radiation_constant", phi);
    rebx_set_param_double(rebx, &yorp->ap, "yorp_lstar", lstar);

    // set parameters for the asteroid
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_rotation_frequency", rotation_frequency);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_c_body", c_body);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_body_density", density);
    particles[1].r = radius;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_tensile_strength", sigma);

    // set obliquity parameters for the asteroid
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_obliquity", obliquity);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_alpha", alpha);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "yorp_beta", beta);

    rebx_add_operator(rebx, yorp); // add YORP effect to the sim

    // initialize output files for plotting
    FILE *fptr1; // time output file
    fptr1 = fopen("times.txt", "w");
    fclose(fptr1);

    FILE *fptr2; // rotation frequency output file
    fptr2 = fopen("omegas.txt", "w");
    fclose(fptr2);

    FILE *fptr3; // obliquity output file
    fptr3 = fopen("obliquities.txt", "w");
    fclose(fptr3);

    // integrate the simulation over tmax time (in this case, 1 Myr)
    double tmax = 1.5E7;
    reb_simulation_integrate(sim, tmax);

    // print final rotation frequency and radius for each asteroid
    for(int i=1; i<(sim->N); i++){
        struct reb_particle* p = &sim->particles[i];

        double* final_w = rebx_get_param(sim->extras, p->ap, "yorp_rotation_frequency"); // final rotation frequency
        *final_w = *final_w/yr_conv; // convert to rad/s

        double* final_eps = rebx_get_param(sim->extras, p->ap, "yorp_obliquity"); // final obliquity
        *final_eps = 180.0*(*final_eps)/M_PI; // convert to degrees

        // print final values
        printf("\nASTEROID %d FINAL ROTATION FREQUENCY: %1.10f rad/s\n", i, *final_w);
        printf("ASTEROID %d FINAL OBLIQUITY: %1.2f deg\n", i, *final_eps);
        printf("ASTEROID %d FINAL RADIUS: %1.5f m\n", i, (p->r)*au_conv);
    }
    rebx_free(rebx); // free REBOUNDx memory
}

// heartbeat function to check for fragmentation due to YORP spin-up and to record data
void heartbeat(struct reb_simulation* sim) {
    tsteps += 1; // increment number of timesteps
    double G = sim->G; // gravitational constant

    // unit conversions
    double au_conv = 1.495978707e11; // meters in one AU
    double yr_conv = 31557600.0; // seconds in one year

    // minimum radius at which a particle will fragment
    double min_r = 10.0/au_conv;

    for(int i=1; i<sim->N; i++){
        struct reb_particle* p = &sim->particles[i];

        // yorp parameters needed to calculate failure spin rate for fragmentation
        const double* density = rebx_get_param(sim->extras, p->ap, "yorp_body_density");
        const double* sigma = rebx_get_param(sim->extras, p->ap, "yorp_tensile_strength");
        double* rotation_frequency = rebx_get_param(sim->extras, p->ap, "yorp_rotation_frequency");
        double r = p->r;

        // record time, rotation frequency, and obliquity of first asteroid every 100 timesteps
        if ((i == 1) && (tsteps % 100 == 0)){
            const double* omega1 = rebx_get_param(sim->extras, p->ap, "yorp_rotation_frequency");
            const double* obliq1 = rebx_get_param(sim->extras, p->ap, "yorp_obliquity");

            FILE *fptr1; // time output file
            fptr1 = fopen("times.txt", "a");
            fprintf(fptr1, "%1.17f\n", sim->t);
            fclose(fptr1);

            FILE *fptr2; // rotation frequency output file
            fptr2 = fopen("omegas.txt", "a");
            fprintf(fptr2, "%1.17f\n", *omega1/yr_conv);
            fclose(fptr2);

            FILE *fptr3; // obliquity output file
            fptr3 = fopen("obliquities.txt", "a");
            fprintf(fptr3, "%1.17f\n", *obliq1);
            fclose(fptr3);
        }

        // Eq. 2 in Veras and Scheeres (2020). Actually the failure spin rate squared.
        double failure_spin_rate = (4.*M_PI*G*(*density))/3. + (2.*(*sigma))/((*density)*r*r)*(2./3.);

        // break apart only if spinning at or above spin rate AND the radius is greater than min fragment size
        if (((*rotation_frequency)*(*rotation_frequency)) >= failure_spin_rate && p->r > min_r){

            // parameters for fragmentation
            double reset_rotation_frequency = 0.0; // reset fragment rotation frequency to 0 after fragmentation
            double primary_mass_ratio = 1.5; // mass ratio of the primary binary component after fragmentation
            double primary_secondary_mass_ratio = 50.0; // mass ratio of the primary binary to the shards created after fragmentation
            int n_shards = 3; // number of shards (small equal-mass fragments) to create

            // reset the rotation frequency
            printf("Asteroid %d break-up at t = %1.5f yr. Rotation frequency at break-up: %1.10f rad/s\n", i, sim->t, *rotation_frequency/yr_conv);
            *rotation_frequency = reset_rotation_frequency;

            // save original particle mass and radius
            double m_old = p->m;
            double r_old = p->r;

            // set mass and radius of first primary fragment
            double m_primary = m_old/(1.0+1.0/(primary_secondary_mass_ratio));
            double r_primary = pow(1.0+1.0/(primary_secondary_mass_ratio), -1.0/3.0)*r_old;
            p->m = m_primary/(1.0+1.0/(primary_mass_ratio));
            p->r = pow(1.0+1.0/(primary_mass_ratio), -1.0/3.0)*r_primary;

            // create second primary fragment and add to simulation
            struct reb_particle fragment = {0};
            fragment.m = m_primary - p->m;
            fragment.r = pow(primary_mass_ratio, -1.0/3.0)*(p->r);
            fragment.x = p->x + 10*(p->r); // offset slightly so not overlapping
            fragment.y = p->y;
            fragment.z = p->z;
            fragment.vx = p->vx;
            fragment.vy = p->vy;
            fragment.vz = p->vz;
            reb_simulation_add(sim, fragment);

            // get other original particle parameters to pass to fragments
            const double* c_body = rebx_get_param(sim->extras, p->ap, "yorp_c_body");
            const double* obliquity = rebx_get_param(sim->extras, p->ap, "yorp_obliquity");
            const double* alpha = rebx_get_param(sim->extras, p->ap, "yorp_alpha");
            const double* beta = rebx_get_param(sim->extras, p->ap, "yorp_beta");

            // address of new particle
            struct reb_particle* new_particle = &sim->particles[sim->N - 1];

            // set simulation-specific parameters for the new fragment
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_rotation_frequency", reset_rotation_frequency);
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_c_body", *c_body);
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_body_density", *density);
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_tensile_strength", *sigma);
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_obliquity", *obliquity);
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_alpha", *alpha);
            rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_beta", *beta);

            // set mass and radius totals for shards
            double m_secondary = m_old - m_primary;
            double r_secondary = r_primary*pow(primary_secondary_mass_ratio, -1.0/3.0);

            // add shards to the simulation
            for(int j=0; j<(n_shards); j++){
                struct reb_particle shard = {0};
                // shards all have the same mass and radius
                shard.m = m_secondary/(n_shards);
                shard.r = pow(n_shards, -1.0/3.0)*r_secondary;
                shard.x = p->x + 10*(j+2)*(p->r); // offset each fragment slightly so they are not overlapping
                shard.y = p->y;
                shard.z = p->z;
                shard.vx = p->vx;
                shard.vy = p->vy;
                shard.vz = p->vz;
                reb_simulation_add(sim, shard);

                // address of new particle
                struct reb_particle* new_particle = &sim->particles[sim->N - 1];

                // set simulation-specific parameters for the new shard
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_rotation_frequency", reset_rotation_frequency);
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_c_body", *c_body);
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_body_density", *density);
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_tensile_strength", *sigma);
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_obliquity", *obliquity);
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_alpha", *alpha);
                rebx_set_param_double(sim->extras, &new_particle->ap, "yorp_beta", *beta);
            }
        }
    }
}