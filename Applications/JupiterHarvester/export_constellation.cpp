///////////////////////////////////////////////////////////////////////////////
// Export Bucket Constellation trajectories as JSON for Three.js viewer
//
// Outputs:
//   constellation_data.json — planet orbits + bucket trajectories + stats
//   Sampled daily for smooth animation
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "jpl_ephemeris.hpp"
#include "lambert.hpp"
#include "../planetary_models/planetary_environment.hpp"

double vec_mag(const double v[3]) { return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
void vec_sub(double o[3], const double a[3], const double b[3]) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
double jd_from_year(double y) { return 2451545.0+(y-2000.0)*365.25; }
double year_from_jd(double jd) { return 2000.0+(jd-2451545.0)/365.25; }
const double AU = 149597870.7; // km

struct LR { double v_inf_d, v_inf_a, tof; double jd_d, jd_a; bool ok; };
LR qlambert(JplEphemeris& eph, JplBody b1, JplBody b2, bool e1, bool e2,
            double jd_lo, double jd_hi, double tof_lo, double tof_hi, double step) {
    LR best; best.ok=false; best.v_inf_d=1e10;
    double sp[3],sv[3];
    for (double jd=jd_lo;jd<=jd_hi;jd+=step) {
        for (double tof=tof_lo;tof<=tof_hi;tof+=step) {
            double p1[3],v1[3],p2[3],v2[3],r1[3],rv1[3],r2[3],rv2[3];
            if(e1) eph.get_earth(jd,p1,v1); else eph.get_state(b1,jd,p1,v1);
            if(e2) eph.get_earth(jd+tof,p2,v2); else eph.get_state(b2,jd+tof,p2,v2);
            eph.get_state(JPL_SUN,jd,sp,sv);
            vec_sub(r1,p1,sp); vec_sub(rv1,v1,sv);
            eph.get_state(JPL_SUN,jd+tof,sp,sv);
            vec_sub(r2,p2,sp); vec_sub(rv2,v2,sv);
            double vt1[3],vt2[3];
            if(!lambert_solve(r1,r2,tof*86400,1.32712440018e11,true,vt1,vt2)) continue;
            double vd[3],va[3]; vec_sub(vd,vt1,rv1); vec_sub(va,vt2,rv2);
            double vid=vec_mag(vd), via=vec_mag(va);
            if(vid<best.v_inf_d && vid<15 && via<15) {
                best.jd_d=jd; best.jd_a=jd+tof; best.v_inf_d=vid;
                best.v_inf_a=via; best.tof=tof; best.ok=true;
            }
        }
    }
    return best;
}

// Bucket state
enum Phase { PH_OUT_V=0, PH_OUT_E=1, PH_OUT_J=2, PH_SCOOP=3, PH_RETURN=4, PH_DONE=5, PH_LOST=6 };
const char* phase_names[] = {"outbound_venus","outbound_earth","outbound_jupiter","scooping","returning","delivered","lost"};

struct BucketTraj {
    int id;
    double jd_launch, jd_venus, jd_earth_flyby, jd_jupiter, jd_depart, jd_home;
    bool valid;
};

int main() {
    fprintf(stderr, "Loading ephemeris...\n");
    JplEphemeris eph;
    if (!eph.load("ephemeris/linux_p1550p2650.440"))
        eph.load("/data/src/zipfel/solar_system/ephemeris/linux_p1550p2650.440");
    if (!eph.is_loaded()) { fprintf(stderr,"FATAL\n"); return 1; }

    double sim_start = 2029.0, sim_end = 2055.0;
    double jd_start = jd_from_year(sim_start), jd_end = jd_from_year(sim_end);
    int BUCKETS_PER_BATCH = 10;
    double SYNODIC = 398.88;

    // Compute all launch window trajectories
    fprintf(stderr, "Computing launch windows...\n");
    std::vector<BucketTraj> templates;
    for (double jd = jd_start; jd < jd_end - 3000; jd += SYNODIC) {
        BucketTraj bt; bt.jd_launch = jd; bt.valid = false;
        LR l1 = qlambert(eph,JPL_EMB,JPL_VENUS,true,false,jd,jd+365,90,250,15);
        if (!l1.ok) { templates.push_back(bt); continue; }
        bt.jd_venus = l1.jd_a;
        LR l2 = qlambert(eph,JPL_VENUS,JPL_EMB,false,true,l1.jd_a,l1.jd_a+500,100,400,15);
        if (!l2.ok) { templates.push_back(bt); continue; }
        bt.jd_earth_flyby = l2.jd_a;
        LR l3 = qlambert(eph,JPL_EMB,JPL_JUPITER,true,false,l2.jd_a,l2.jd_a+400,500,1200,30);
        if (!l3.ok) { templates.push_back(bt); continue; }
        bt.jd_jupiter = l3.jd_a;
        double jd_filled = l3.jd_a + 620; // scoop campaign
        LR l4 = qlambert(eph,JPL_JUPITER,JPL_EMB,false,true,jd_filled,jd_filled+800,500,1200,30);
        if (!l4.ok) { templates.push_back(bt); continue; }
        bt.jd_depart = l4.jd_d;
        bt.jd_home = l4.jd_a;
        bt.valid = true;
        templates.push_back(bt);
        fprintf(stderr, "  Window %.1f: valid, loop=%.1f yr\n",
                year_from_jd(jd), (bt.jd_home-jd)/365.25);
    }

    // Create bucket fleet
    struct Bucket {
        int id; BucketTraj traj; bool alive; double he_delivered;
    };
    std::vector<Bucket> fleet;
    int next_id = 1;
    for (auto& bt : templates) {
        if (!bt.valid) continue;
        for (int b = 0; b < BUCKETS_PER_BATCH; b++) {
            Bucket bk; bk.id = next_id++; bk.traj = bt;
            bk.alive = ((rand()%100) >= 10); // 10% attrition
            bk.he_delivered = bk.alive ? 8250 : 0; // kg He per successful loop
            fleet.push_back(bk);
        }
    }

    fprintf(stderr, "Fleet: %d buckets\n", (int)fleet.size());

    // Output JSON
    FILE* fp = fopen("constellation_data.json", "w");
    if (!fp) { fprintf(stderr, "Cannot open output\n"); return 1; }

    fprintf(fp, "{\n");

    // Metadata
    fprintf(fp, "  \"title\": \"THE BUCKET — Jupiter Helium Constellation\",\n");
    fprintf(fp, "  \"start_year\": %.1f,\n", sim_start);
    fprintf(fp, "  \"end_year\": %.1f,\n", sim_end);
    fprintf(fp, "  \"au_km\": %.1f,\n", AU);
    fprintf(fp, "  \"bucket_count\": %d,\n", (int)fleet.size());

    // Planet orbits (sampled every 10 days)
    fprintf(fp, "  \"planets\": {\n");
    const char* pnames[] = {"venus","earth","mars","jupiter"};
    JplBody pbodies[] = {JPL_VENUS, JPL_EMB, JPL_MARS, JPL_JUPITER};
    bool pearth[] = {false, true, false, false};

    for (int p = 0; p < 4; p++) {
        fprintf(fp, "    \"%s\": [\n", pnames[p]);
        bool first = true;
        for (double jd = jd_start; jd <= jd_end; jd += 10) {
            double pos[3], vel[3], sp[3], sv[3];
            if (pearth[p]) eph.get_earth(jd, pos, vel);
            else eph.get_state(pbodies[p], jd, pos, vel);
            eph.get_state(JPL_SUN, jd, sp, sv);
            double x = (pos[0]-sp[0])/AU, y = (pos[1]-sp[1])/AU, z = (pos[2]-sp[2])/AU;
            if (!first) fprintf(fp, ",\n");
            fprintf(fp, "      [%.4f, %.6f, %.6f, %.6f]", year_from_jd(jd), x, y, z);
            first = false;
        }
        fprintf(fp, "\n    ]%s\n", p < 3 ? "," : "");
    }
    fprintf(fp, "  },\n");

    // Bucket trajectories
    fprintf(fp, "  \"buckets\": [\n");
    for (int bi = 0; bi < (int)fleet.size(); bi++) {
        Bucket& bk = fleet[bi];
        BucketTraj& t = bk.traj;
        fprintf(fp, "    {\n");
        fprintf(fp, "      \"id\": %d,\n", bk.id);
        fprintf(fp, "      \"alive\": %s,\n", bk.alive ? "true" : "false");
        fprintf(fp, "      \"he_delivered\": %.0f,\n", bk.he_delivered);

        // Keyframes: position at each phase transition
        // Interpolate position along transfer arcs using planet positions
        fprintf(fp, "      \"waypoints\": [\n");

        // Helper: get heliocentric position of a body at a time
        auto get_pos = [&](JplBody body, bool is_earth, double jd, double& ox, double& oy, double& oz) {
            double pos[3], vel[3], sp[3], sv[3];
            if (is_earth) eph.get_earth(jd, pos, vel);
            else eph.get_state(body, jd, pos, vel);
            eph.get_state(JPL_SUN, jd, sp, sv);
            ox = (pos[0]-sp[0])/AU; oy = (pos[1]-sp[1])/AU; oz = (pos[2]-sp[2])/AU;
        };

        struct WP { double year, x, y, z; const char* phase; double tank_pct; };
        std::vector<WP> wps;

        double ex,ey,ez, vx,vy,vz, jx,jy,jz;

        // Launch (at Earth)
        get_pos(JPL_EMB, true, t.jd_launch, ex,ey,ez);
        wps.push_back({year_from_jd(t.jd_launch), ex,ey,ez, "launch", 0});

        // Venus flyby
        get_pos(JPL_VENUS, false, t.jd_venus, vx,vy,vz);
        wps.push_back({year_from_jd(t.jd_venus), vx,vy,vz, "venus_flyby", 0});

        // Earth flyby
        get_pos(JPL_EMB, true, t.jd_earth_flyby, ex,ey,ez);
        wps.push_back({year_from_jd(t.jd_earth_flyby), ex,ey,ez, "earth_flyby", 0});

        // Jupiter arrival
        get_pos(JPL_JUPITER, false, t.jd_jupiter, jx,jy,jz);
        wps.push_back({year_from_jd(t.jd_jupiter), jx,jy,jz, "jupiter_arrive", 0});

        // Scooping (stays at Jupiter, tank fills)
        for (double f = 0.25; f <= 1.0; f += 0.25) {
            double jd_s = t.jd_jupiter + f * 620;
            get_pos(JPL_JUPITER, false, jd_s, jx,jy,jz);
            // Tiny offset so buckets aren't all on top of Jupiter
            double offset = (bk.id % 20 - 10) * 0.002;
            wps.push_back({year_from_jd(jd_s), jx+offset,jy+offset,jz,
                          "scooping", f*100});
        }

        // Jupiter departure
        get_pos(JPL_JUPITER, false, t.jd_depart, jx,jy,jz);
        wps.push_back({year_from_jd(t.jd_depart), jx,jy,jz, "depart_jupiter", 100});

        // Return: midpoint (interpolate between Jupiter and Earth positions)
        double jd_mid = (t.jd_depart + t.jd_home) / 2;
        get_pos(JPL_JUPITER, false, t.jd_depart, jx,jy,jz);
        get_pos(JPL_EMB, true, t.jd_home, ex,ey,ez);
        wps.push_back({year_from_jd(jd_mid), (jx+ex)/2,(jy+ey)/2,(jz+ez)/2, "returning", 100});

        // Earth arrival
        get_pos(JPL_EMB, true, t.jd_home, ex,ey,ez);
        if (bk.alive)
            wps.push_back({year_from_jd(t.jd_home), ex,ey,ez, "delivered", 100});
        else
            wps.push_back({year_from_jd(t.jd_home), ex,ey,ez, "lost", 0});

        for (int wi = 0; wi < (int)wps.size(); wi++) {
            WP& w = wps[wi];
            fprintf(fp, "        {\"year\": %.4f, \"x\": %.6f, \"y\": %.6f, \"z\": %.6f, \"phase\": \"%s\", \"tank\": %.0f}%s\n",
                    w.year, w.x, w.y, w.z, w.phase, w.tank_pct,
                    wi < (int)wps.size()-1 ? "," : "");
        }

        fprintf(fp, "      ]\n");
        fprintf(fp, "    }%s\n", bi < (int)fleet.size()-1 ? "," : "");
    }
    fprintf(fp, "  ],\n");

    // Cumulative stats by year
    fprintf(fp, "  \"stats\": [\n");
    double cum_he = 0;
    int cum_delivered = 0, cum_lost = 0, cum_launched = 0;
    bool first_stat = true;
    for (double yr = sim_start; yr <= sim_end; yr += 0.5) {
        int launched = 0, delivered = 0, lost = 0, active = 0;
        double he_yr = 0;
        for (auto& bk : fleet) {
            if (year_from_jd(bk.traj.jd_launch) <= yr) {
                launched++;
                if (year_from_jd(bk.traj.jd_home) <= yr) {
                    if (bk.alive) { delivered++; he_yr += bk.he_delivered; }
                    else lost++;
                } else {
                    active++;
                }
            }
        }
        cum_launched = launched;
        cum_delivered = delivered;
        cum_lost = lost;
        cum_he = he_yr;

        if (!first_stat) fprintf(fp, ",\n");
        fprintf(fp, "    {\"year\": %.1f, \"launched\": %d, \"active\": %d, \"delivered\": %d, \"lost\": %d, \"he_tonnes\": %.1f}",
                yr, launched, active, delivered, lost, cum_he/1000);
        first_stat = false;
    }
    fprintf(fp, "\n  ]\n");

    fprintf(fp, "}\n");
    fclose(fp);

    fprintf(stderr, "Wrote constellation_data.json\n");
    fprintf(stderr, "Buckets: %d, Valid trajectories: %d\n",
            (int)fleet.size(), (int)templates.size());
    return 0;
}
