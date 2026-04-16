#ifndef CATALOGS_H
#define CATALOGS_H

// constants

#define CONFIG_SIZE 4

#define N_STARS 1606
#define N_PAIRS 330309
#define HASH_TABLE_SIZE 1200
#define MAX_ANGLE 60
#define DELTA_HASH 0.05

#define PI 3.14159265

// settings

#define DEG_PER_PIXEL 0.027778 // 0.027766
#define DIST_ERROR 2
//#define MAX_SQ_ERR (DIST_ERROR * DIST_ERROR * CONFIG_SIZE)
#define MAX_SQ_ERR 0.1

// star detection parameters

#define MIN_STAR_SIZE 5
#define MAX_STAR_SIZE 300
#define MIN_PIXEL_BRIGHTNESS 30
#define MIN_STAR_BRIGHTNESS 500
#define MAX_STAR_BRIGHTNESS 80000



extern const double const RA[N_STARS];

extern const double const dec[N_STARS];

extern const double const dist_mat[N_STARS * (N_STARS - 1) / 2];

extern const int const k[N_PAIRS];

extern const int const l[N_PAIRS];

extern const double const d[N_PAIRS];

extern const int const hash_table[HASH_TABLE_SIZE];




#endif