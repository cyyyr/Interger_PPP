//
// Created by cyr on 03.11.2020.
//


#include "Positioning.h"

/* single to double-difference transformation matrix (D') --------------------*/
int Positioning::ddmat(data_t data, Matrix<double> D) {
    int i,j,k,m,f;
    int nb = 0;
    int nx = data.nx;
    int na = data.na;
    int nf = NF_RTK(&data.opt);
    int nofix;

    for (i = 0; i < MAXSAT; i++) {
        for (j = 0; j < NFREQ; j++) {
            data->ssat[i].fix[j] = 0;
        }
    }
    for (i = 0; i < na; i++) {
        D(i, i) = 1.0;
    }

    for (m = 0; m < 4; m++) { /* m=0:gps/qzs/sbs, 1:glo, 2:gal, 3:bds */
        nofix = (m == 1 && rtk->opt.glomodear == 0) || (m == 3 && rtk->opt.bdsmodear == 0);

        for (f = 0, k = na; f < nf; f++, k += MAXSAT) {
            if (i < k + MAXSAT) {
                for (i = k; i < k + MAXSAT; i++) {
                    if (rtk->x[i] == 0.0 || !test_sys(rtk->ssat[i - k].sys, m) ||
                        !rtk->ssat[i - k].vsat[f] || !rtk->ssat[i - k].half[f]) {
                        continue;
                    }
                    if (rtk->ssat[i - k].lock[f] > 0 && !(rtk->ssat[i - k].slip[f] & 2) &&
                        rtk->ssat[i - k].azel[1] >= rtk->opt.elmaskar && !nofix) {
                        rtk->ssat[i - k].fix[f] = 2; /* fix */
                        break;
                    }
                    rtk->ssat[i - k].fix[f] = 1;
                }
                for (j = k; j < k + MAXSAT; j++) {
                    if (i == j || rtk->x[j] == 0.0 || !test_sys(rtk->ssat[j - k].sys, m) ||
                        !rtk->ssat[j - k].vsat[f]) {
                        continue;
                    }
                    if (rtk->ssat[j - k].lock[f] > 0 && !(rtk->ssat[j - k].slip[f] & 2) &&
                        rtk->ssat[i - k].vsat[f] &&
                        rtk->ssat[j - k].azel[1] >= rtk->opt.elmaskar && !nofix) {
                        D[i + (na + nb) * nx] = 1.0;
                        D[j + (na + nb) * nx] = -1.0;
                        nb++;
                        rtk->ssat[j - k].fix[f] = 2; /* fix */
                    } else {
                        rtk->ssat[j - k].fix[f] = 1;
                    }
                }
            }
        }
    }
    trace(5, "D=\n");
    tracemat(5, D, nx, na + nb, 2, 0);
    return nb;
}

/* restore single-differenced ambiguity --------------------------------------*/
void restamb(rtk_t *rtk, const double *bias, int nb __attribute((unused)), double *xa) {
    int i;
    int n;
    int m;
    int f;
    int index[MAXSAT];
    int nv = 0;
    int nf = NF_RTK(&rtk->opt);

    trace(3, "restamb :\n");

    for (i = 0; i < rtk->nx; i++) {
        xa[i] = rtk->x[i];
    }
    for (i = 0; i < rtk->na; i++) {
        xa[i] = rtk->xa[i];
    }

    for (m = 0; m < 4; m++) {
        for (f = 0; f < nf; f++) {
            for (n = i = 0; i < MAXSAT; i++) {
                if (!test_sys(rtk->ssat[i].sys, m) || rtk->ssat[i].fix[f] != 2) {
                    continue;
                }
                index[n++] = IB_RTK(i + 1, f, &rtk->opt);
            }
            if (n < 2) {
                continue;
            }

            xa[index[0]] = rtk->x[index[0]];

            for (i = 1; i < n; i++) {
                xa[index[i]] = xa[index[0]] - bias[nv++];
            }
        }
    }
}


/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
int resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa) {
    prcopt_t *opt = &rtk->opt;
    int i;
    int j;
    int ny;
    int nb;
    int info;
    int nx = rtk->nx;
    int na = rtk->na;
    double *D;
    double *DP;
    double *y;
    double *Qy;
    double *b;
    double *db;
    double *Qb;
    double *Qab;
    double *QQ;
    double s[2];

    trace(3, "resamb_LAMBDA : nx=%d\n", nx);

    rtk->sol.ratio = 0.0;

    if (rtk->opt.mode <= PMODE_DGPS || rtk->opt.modear == ARMODE_OFF ||
        rtk->opt.thresar[0] < 1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D = zeros(nx, nx);
    if ((nb = ddmat(rtk, D)) <= 0) {
        errmsg(rtk, "no valid double-difference\n");
        free(D);
        return 0;
    }
    ny = na + nb;
    y = mat(ny, 1);
    Qy = mat(ny, ny);
    DP = mat(ny, nx);
    b = mat(nb, 2);
    db = mat(nb, 1);
    Qb = mat(nb, nb);
    Qab = mat(na, nb);
    QQ = mat(na, nb);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    matmul("TN", ny, 1, nx, 1.0, D, rtk->x, 0.0, y);
    matmul("TN", ny, nx, nx, 1.0, D, rtk->P, 0.0, DP);
    matmul("NN", ny, ny, nx, 1.0, DP, D, 0.0, Qy);

    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i = 0; i < nb; i++) {
        for (j = 0; j < nb; j++) {
            Qb[i + j * nb] = Qy[na + i + (na + j) * ny];
        }
    }
    for (i = 0; i < na; i++) {
        for (j = 0; j < nb; j++) {
            Qab[i + j * na] = Qy[i + (na + j) * ny];
        }
    }

    trace(4, "N(0)=");
    tracemat(4, y + na, 1, nb, 10, 3);

    /* lambda/mlambda integer least-square estimation */
    if (!(info = lambda(nb, 2, y + na, Qb, b, s))) {
        trace(4, "N(1)=");
        tracemat(4, b, 1, nb, 10, 3);
        trace(4, "N(2)=");
        tracemat(4, b + nb, 1, nb, 10, 3);

        rtk->sol.ratio = s[0] > 0 ? static_cast<float>(s[1] / s[0]) : 0.0F;
        if (rtk->sol.ratio > 999.9) {
            rtk->sol.ratio = 999.9F;
        }

        /* validation by popular ratio-test */
        if (s[0] <= 0.0 || s[1] / s[0] >= opt->thresar[0]) {
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i = 0; i < na; i++) {
                rtk->xa[i] = rtk->x[i];
                for (j = 0; j < na; j++) {
                    rtk->Pa[i + j * na] = rtk->P[i + j * nx];
                }
            }
            for (i = 0; i < nb; i++) {
                bias[i] = b[i];
                y[na + i] -= b[i];
            }
            if (!matinv(Qb, nb)) {
                matmul("NN", nb, 1, nb, 1.0, Qb, y + na, 0.0, db);
                matmul("NN", na, 1, nb, -1.0, Qab, db, 1.0, rtk->xa);

                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN", na, nb, nb, 1.0, Qab, Qb, 0.0, QQ);
                matmul("NT", na, na, nb, -1.0, QQ, Qab, 1.0, rtk->Pa);

                trace(3, "resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb, s[0] == 0.0 ? 0.0 : s[1] / s[0], s[0], s[1]);

                /* restore single-differenced ambiguity */
                restamb(rtk, bias, nb, xa);
            } else {
                nb = 0;
            }
        } else { /* validation failed */
            errmsg(rtk, "ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb, s[1] / s[0], s[0], s[1]);
            nb = 0;
        }
    } else {
        errmsg(rtk, "lambda error (info=%d)\n", info);
    }
    free(D);
    free(y);
    free(Qy);
    free(DP);
    free(b);
    free(db);
    free(Qb);
    free(Qab);
    free(QQ);

    return nb; /* number of ambiguities */
}

/* validation of solution ----------------------------------------------------*/
int valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
           int nv, double thres) {
#if 0
    prcopt_t *opt = &rtk->opt;
    double vv = 0.0;
#endif
    double fact = thres * thres;
    int i;
    int stat = 1;
    int sat1;
    int sat2;
    int type;
    int freq;
    char stype;

    trace(3, "valpos  : nv=%d thres=%.1f\n", nv, thres);

    /* post-fit residual test */
    for (i = 0; i < nv; i++) {
        if (v[i] * v[i] <= fact * R[i + i * nv]) {
            continue;
        }
        sat1 = (vflg[i] >> 16) & 0xFF;
        sat2 = (vflg[i] >> 8) & 0xFF;
        type = (vflg[i] >> 4) & 0xF;
        freq = vflg[i] & 0xF;
        stype = type == 0 ? 'L' : (type == 1 ? 'L' : 'C');
        errmsg(rtk, "large residual (sat=%2d-%2d %c%d v=%6.3f sig=%.3f)\n",
               sat1, sat2, stype, freq + 1, v[i], std::sqrt(R[i + i * nv]));
    }
#if 0 /* omitted v.2.4.0 */
    if (stat && nv > NP(opt))
        {
            /* chi-square validation */
            for (i = 0; i < nv; i++) vv += v[i] * v[i] / R[i + i * nv];

            if (vv > chisqr[nv - NP(opt) - 1])
                {
                    errmsg(rtk, "residuals validation failed (nv=%d np=%d vv=%.2f cs=%.2f)\n",
                        nv, NP(opt), vv, chisqr[nv - NP(opt) - 1]);
                    stat = 0;
                }
            else
                {
                    trace(3, "valpos : validation ok (%s nv=%d np=%d vv=%.2f cs=%.2f)\n",
                        rtk->tstr, nv, NP(opt), vv, chisqr[nv - NP(opt) - 1]);
                }
        }
#endif
    return stat;
}