// Copyright (c) 2026 Darrell Thomas. MIT License.
// See LICENSE file for details.

///////////////////////////////////////////////////////////////////////////////
// JPL Development Ephemeris (DE440) Reader
//
// Reads JPL Linux binary ephemeris files and returns planet positions/velocities
// at any Julian Date within the file's coverage.
//
// DE440 Linux binary format (rec_size = 8144 bytes = 1018 doubles):
//   Record 1 (bytes 0-8143): Title + constant names + pointers
//     0-251: 3×84-char title lines
//     252-2651: 400 constant names (6 chars each)
//     2652-2675: JD_start, JD_end, block_interval (3 doubles)
//     2676-2679: NCON (int)
//     2680-2695: AU, EMRAT (2 doubles)
//     2696-2839: 12 pointer triplets (12×3 ints = 144 bytes)
//     2840-2843: DE number (int)
//     2844-2855: 13th pointer triplet (librations)
//   Record 2 (bytes 8144-16287): Constant values (NCON doubles)
//   Record 3+ (bytes 16288+): Data blocks (each 1018 doubles)
//     First 2 doubles per block: JD_start, JD_end of block
//     Remaining: Chebyshev coefficients per pointer triplets
//
// Reference: Park et al. 2021, AJ 161:105
// 20260413 Created for Jupiter Harvester / CADAC++ solar system propagator
///////////////////////////////////////////////////////////////////////////////

#ifndef JPL_EPHEMERIS_HPP
#define JPL_EPHEMERIS_HPP

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

enum JplBody {
    JPL_MERCURY = 1, JPL_VENUS = 2, JPL_EMB = 3, JPL_MARS = 4,
    JPL_JUPITER = 5, JPL_SATURN = 6, JPL_URANUS = 7, JPL_NEPTUNE = 8,
    JPL_PLUTO = 9, JPL_MOON = 10, JPL_SUN = 11,
    JPL_NUTATION = 12, JPL_LIBRATION = 13
};

class JplEphemeris
{
public:
    JplEphemeris() : fp_(nullptr), loaded_(false), rec_size_(8144) {}
    ~JplEphemeris() { close(); }

    bool load(const char* filepath)
    {
        close();
        fp_ = fopen(filepath, "rb");
        if (!fp_) {
            fprintf(stderr, "JplEphemeris: Cannot open %s\n", filepath);
            return false;
        }

        // Record 1: title + names + pointers
        // JD range at offset 2652
        fseek(fp_, 2652, SEEK_SET);
        fread(&jd_start_, sizeof(double), 1, fp_);
        fread(&jd_end_, sizeof(double), 1, fp_);
        fread(&block_interval_, sizeof(double), 1, fp_);

        // NCON at offset 2676
        fread(&ncon_, sizeof(int), 1, fp_);

        // AU and EMRAT at offset 2680
        fread(&au_, sizeof(double), 1, fp_);
        fread(&emrat_, sizeof(double), 1, fp_);

        // 12 pointer triplets at offset 2696
        for (int i = 0; i < 12; i++) {
            fread(&ipt_[i][0], sizeof(int), 1, fp_);
            fread(&ipt_[i][1], sizeof(int), 1, fp_);
            fread(&ipt_[i][2], sizeof(int), 1, fp_);
        }

        // DE number at offset 2840
        fread(&de_num_, sizeof(int), 1, fp_);

        // 13th pointer triplet (librations) at offset 2844
        fread(&ipt_[12][0], sizeof(int), 1, fp_);
        fread(&ipt_[12][1], sizeof(int), 1, fp_);
        fread(&ipt_[12][2], sizeof(int), 1, fp_);

        // Record size = 1018 doubles = 8144 bytes (DE440)
        ncoeff_ = rec_size_ / (int)sizeof(double); // 1018
        data_offset_ = 2 * rec_size_; // data starts at record 3

        loaded_ = true;
        printf("JplEphemeris: Loaded DE%d — JD %.1f to %.1f (%.0f to %.0f CE), "
               "interval=%.0f days, AU=%.3f km, EMRAT=%.6f\n",
               de_num_, jd_start_, jd_end_,
               2000.0 + (jd_start_ - 2451545.0)/365.25,
               2000.0 + (jd_end_ - 2451545.0)/365.25,
               block_interval_, au_, emrat_);

        return true;
    }

    void close() { if (fp_) { fclose(fp_); fp_ = nullptr; } loaded_ = false; }
    bool is_loaded() const { return loaded_; }
    double jd_start() const { return jd_start_; }
    double jd_end() const { return jd_end_; }
    int de_number() const { return de_num_; }
    double au_km() const { return au_; }
    double emrat() const { return emrat_; }

    // Get position (km) and velocity (km/s) of a body at Julian Date
    // All coordinates are Solar System Barycentric (SSB), J2000
    bool get_state(JplBody body, double jd, double pos[3], double vel[3]) const
    {
        if (!loaded_ || jd < jd_start_ || jd > jd_end_) return false;

        int idx = static_cast<int>(body) - 1;
        if (idx < 0 || idx > 12) return false;

        double pv[6];

        if (body == JPL_MOON) {
            // Moon is stored geocentric; convert to barycentric
            // Moon_bary = EMB + Moon_geo * EMRAT/(1+EMRAT)
            double emb_pv[6], mgeo_pv[6];
            interpolate(2, jd, emb_pv);  // EMB (index 2)
            interpolate(9, jd, mgeo_pv); // Moon geocentric (index 9)
            double fac = emrat_ / (1.0 + emrat_);
            for (int i = 0; i < 6; i++)
                pv[i] = emb_pv[i] + mgeo_pv[i] * fac;
        } else {
            interpolate(idx, jd, pv);
        }

        pos[0] = pv[0]; pos[1] = pv[1]; pos[2] = pv[2];
        // Velocity: interpolate returns km/day, convert to km/s
        vel[0] = pv[3] / 86400.0;
        vel[1] = pv[4] / 86400.0;
        vel[2] = pv[5] / 86400.0;
        return true;
    }

    // Get Earth position (derived from EMB - Moon/(1+EMRAT))
    bool get_earth(double jd, double pos[3], double vel[3]) const
    {
        double emb_pv[6], mgeo_pv[6];
        interpolate(2, jd, emb_pv);  // EMB
        interpolate(9, jd, mgeo_pv); // Moon geocentric
        double fac = 1.0 / (1.0 + emrat_);
        pos[0] = emb_pv[0] - mgeo_pv[0] * fac;
        pos[1] = emb_pv[1] - mgeo_pv[1] * fac;
        pos[2] = emb_pv[2] - mgeo_pv[2] * fac;
        vel[0] = (emb_pv[3] - mgeo_pv[3] * fac) / 86400.0;
        vel[1] = (emb_pv[4] - mgeo_pv[4] * fac) / 86400.0;
        vel[2] = (emb_pv[5] - mgeo_pv[5] * fac) / 86400.0;
        return true;
    }

private:
    FILE* fp_;
    bool loaded_;
    int rec_size_;
    double jd_start_, jd_end_, block_interval_;
    double au_, emrat_;
    int ncon_, de_num_, ncoeff_;
    int ipt_[13][3];
    long data_offset_;

    // Interpolate Chebyshev polynomials for body_idx (0-based)
    // Returns position (km) in pv[0..2] and velocity (km/day) in pv[3..5]
    void interpolate(int body_idx, double jd, double pv[6]) const
    {
        int ncomp = (body_idx == 11) ? 2 : 3; // nutations have 2 components
        int offset = ipt_[body_idx][0] - 1;     // 0-based offset in record (doubles)
        int ncoef  = ipt_[body_idx][1];
        int nsub   = ipt_[body_idx][2];

        if (ncoef == 0 || nsub == 0) {
            for (int i = 0; i < 6; i++) pv[i] = 0;
            return;
        }

        // Which data block contains this JD?
        int block_num = (int)((jd - jd_start_) / block_interval_);
        if (block_num < 0) block_num = 0;
        double block_start = jd_start_ + block_num * block_interval_;

        // Which subinterval within the block?
        double sub_interval = block_interval_ / nsub;
        int sub_num = (int)((jd - block_start) / sub_interval);
        if (sub_num >= nsub) sub_num = nsub - 1;
        if (sub_num < 0) sub_num = 0;
        double sub_start = block_start + sub_num * sub_interval;

        // Normalized time in [-1, 1]
        double tc = 2.0 * (jd - sub_start) / sub_interval - 1.0;
        double vfac = 2.0 / sub_interval; // velocity conversion: d(tc)/d(t) in 1/day

        // File offset: skip 2 header records + block_num data records + 2 JD doubles + coefficients
        // Each data record starts with JD_start, JD_end (2 doubles), then coefficients
        long file_pos = data_offset_ + (long)block_num * rec_size_
                      + (long)(offset + sub_num * ncomp * ncoef) * sizeof(double);

        // Read coefficients for all components
        double coef[3 * 20]; // max 3 components × 20 coefficients
        int total_coef = ncomp * ncoef;
        fseek(fp_, file_pos, SEEK_SET);
        fread(coef, sizeof(double), total_coef, fp_);

        // Evaluate Chebyshev polynomials
        for (int comp = 0; comp < ncomp && comp < 3; comp++)
        {
            double* c = coef + comp * ncoef;

            // T_0 = 1, T_1 = tc
            double T0 = 1.0, T1 = tc;
            double pos_val = c[0] * T0 + ((ncoef > 1) ? c[1] * T1 : 0.0);

            // Derivatives: dT_0/dtc = 0, dT_1/dtc = 1
            double dT0 = 0.0, dT1 = 1.0;
            double vel_val = (ncoef > 1) ? c[1] : 0.0;

            for (int n = 2; n < ncoef; n++)
            {
                double T2 = 2.0 * tc * T1 - T0;
                double dT2 = 2.0 * (T1 + tc * dT1) - dT0;

                pos_val += c[n] * T2;
                vel_val += c[n] * dT2;

                T0 = T1; T1 = T2;
                dT0 = dT1; dT1 = dT2;
            }

            pv[comp] = pos_val;
            pv[comp + 3] = vel_val * vfac; // km/day
        }

        for (int comp = ncomp; comp < 3; comp++) {
            pv[comp] = 0; pv[comp + 3] = 0;
        }
    }
};

#endif // JPL_EPHEMERIS_HPP
