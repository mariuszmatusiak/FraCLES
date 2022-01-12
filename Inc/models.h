/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    models.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   models.h header file with example models obtained in MATLAB  (paper
 *           Matusiak2020l), including: ORA approximation of FOTF, balred
 *          reduced ORA approximation, first order plant, second order plant,
 *          fractional order plant
 * @date    2018/09/24
 ******************************************************************************/

#ifndef MODELS_H_
#define MODELS_H_

#include "types.h"

// Example plants

/**
 * \brief fotf_best
 *
 * Details:
 * Fractional-order transfer function:
 *                      1
 * -------------------------------------------- exp(-1s)
 * 0.18234s^{1.9909}+0.65536s^{0.98319}+0.99929
 */
const fotf fotf_best = {
        .num_order_length = 1,
        .den_order_length = 3,
        .delay = 1
};
const floatp_t fotf_best_b[] =  {1};
const floatp_t fotf_best_nb[] = {0};
const floatp_t fotf_best_a[] =  {0.18234, 0.65536, 0.99929};
const floatp_t fotf_best_na[] = {1.99090, 0.98319, 0};


/**
 * \brief  fotf_best_oustaloup_3_balred_3_discrete
 * Details:
 *
 *             0.00142 - 0.004268 z^-1 + 0.004281 z^-2 - 0.001433 z^-3
 * z^(-1000) * -------------------------------------------------------
 *                    1 - 2.996 z^-1 + 2.992 z^-2 - 0.9959 z^-3
 *
 * Sample time: 0.001 seconds
 * Discrete-time transfer function.
 */
const discretePlantCoeffs fotf_best_oustaloup_3_balred_3_discrete = {
        .b[0] = (floatp_t) 0.00141997809862250,
        .b[1] = (floatp_t)-0.00426750404873265,
        .b[2] = (floatp_t) 0.00428081749129811,
        .b[3] = (floatp_t)-0.00143328952852651,
        .nb[0] = 0,
        .nb[1] = -1,
        .nb[2] = -2,
        .nb[3] = -3,
        .a[0]  = (floatp_t) 1.0,
        .a[1]  = (floatp_t)-2.99589971448314,
        .a[2]  = (floatp_t) 2.99180655152227,
        .a[3]  = (floatp_t)-0.995906835027740,
        .na[0] = 0,
        .na[1] = -1,
        .na[2] = -2,
        .na[3] = -3,
        .num_order = 3,
        .den_order = 3,
        .delay = 0
};

const discretePlantCoeffs dp1d = {
        .b[0] = (floatp_t)0.000123629474891807,
        .b[1] = (floatp_t)0.00190115587019793,
        .nb[0] = 0,
        .nb[1] = -1,
        .a[0]  = 1.f,
        .a[1]  = (floatp_t)-0.997975214654910,
        .na[0] = 0,
        .na[1] = -1,
        .num_order = 1,
        .den_order = 1,
        .delay = 1228
};

const discretePlantCoeffs dp2d = {
        .b[0] = (floatp_t)2.90801788577009e-06,
        .b[1] = (floatp_t)6.60568256821855e-06,
        .b[2] = (floatp_t)2.55473874360670e-07,
        .nb[0] = 0,
        .nb[1] = -1,
        .nb[2] = -2,
        .a[0]  = 1.f,
        .a[1]  = (floatp_t)-1.99374886431811,
        .a[2]  = (floatp_t)0.993758633492438,
        .na[0] = 0,
        .na[1] = -1,
        .na[2] = -2,
        .num_order = 2,
        .den_order = 2,
        .delay = 1065
};

const discretePlantCoeffs dpfd_zoh = {
        .b[0] = (floatp_t)-0.028242534341347499,
        .b[1] = (floatp_t) 0.085449792586368201,
        .b[2] = (floatp_t)-0.086177578523294704,
        .b[3] = (floatp_t) 0.0289703554269399,
        .nb[0] = 0,
        .nb[1] = -1,
        .nb[2] = -2,
        .nb[3] = -3,
        .a[0]  = (floatp_t) 1.0,
        .a[1]  = (floatp_t)-2.9909770609612689,
        .a[2]  = (floatp_t) 2.9819824191408282,
        .a[3]  = (floatp_t)-0.991005323617552,
        .na[0] = 0,
        .na[1] = -1,
        .na[2] = -2,
        .na[3] = -3,
        .num_order = 3,
        .den_order = 3,
        .delay = 700
};

const discretePlantCoeffs dpfd_tustin = {
        .b[0] = (floatp_t)-0.027754019,
        .b[1] =  (floatp_t)0.083981858,
        .b[2] = (floatp_t)-0.084707238,
        .b[3] =  (floatp_t)0.028479434,
        .nb[0] = 0,
        .nb[1] = -1,
        .nb[2] = -2,
        .nb[3] = -3,
        .a[0]  = 1.f,
        .a[1]  =  (floatp_t)2.99097705,
        .a[2]  = (floatp_t)-2.98198241,
        .a[3]  =  (floatp_t)0.99100532,
        .na[0] = 0,
        .na[1] = -1,
        .na[2] = -2,
        .na[3] = -3,
        .num_order = 3,
        .den_order = 3,
        .delay = 0
};


#endif /* MODELS_H_ */
