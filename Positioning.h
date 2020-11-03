//
// Created by cyr on 03.11.2020.
//

#pragma once
#ifndef ILS_PPP_POSITIONING_H
#define ILS_PPP_POSITIONING_H

#include "Matrix.h"
#include "Lambda.h"


class Positioning {
private:

    struct gtime_t { /* time struct */
        time_t time; /* time (s) expressed by standard time_t */
        double sec;  /* fraction of second under 1 s */
    };

    struct ambc_t {         /* ambiguity control type */
        gtime_t epoch[4];   /* last epoch */
        int n[4];           /* number of epochs */
        double LC[4];       /* linear combination average */
        double LCv[4];      /* linear combination variance */
        int fixcnt;         /* fix count */
    };

    struct sol_t
    {                       /* solution type */
        gtime_t time;       /* time (GPST) */
        double rr[6];       /* position/velocity (m|m/s) */
        /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
        float qr[6];        /* position variance/covariance (m^2) */
        /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
        /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
        double dtr[6];      /* receiver clock bias to time systems (s) */
        unsigned char type; /* type (0:xyz-ecef,1:enu-baseline) */
        unsigned char stat; /* solution status (SOLQ_???) */
        unsigned char ns;   /* number of valid satellites */
        float age;          /* age of differential (s) */
        float ratio;        /* AR ratio factor for validation */
        float thres;        /* AR ratio threshold for validation */
    };

    struct data_t {             /* RTK control/result type */
        sol_t sol;              /* RTK solution */
        double rb[6];           /* base position/velocity (ecef) (m|m/s) */
        int nx, na;             /* number of float states/fixed states */
        double tt;              /* time difference between current and previous (s) */
        double *x, *P;          /* float states and their covariance */
        double *xa, *Pa;        /* fixed states and their covariance */
        int nfix;               /* number of continuous fixes of ambiguity */
        std::vector<ambc_t> ambc;    /* ambiguity control */
    };

    int ddmat(std::vector<data_t> data, Matrix<double> D);

};


#endif //ILS_PPP_POSITIONING_H
