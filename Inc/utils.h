/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    utils.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   utils.h header file with useful functions declarations for FraCLES
 * @date    2019/07/07
 ******************************************************************************/

#ifndef INC_UTILS_H_
#define INC_UTILS_H_

#include "types.h"
#include "../Lib/kiss_fft130/kiss_fft.h"

/**
 *\def LIMIT
 */
#ifndef LIMIT
#define LIMIT(x, min, max)  x = ((x)>(max))?(max):((x) < (min))?(min):(x)
#endif
/**
 * \def MAX
 */
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
/**
 * \def MIN
 */
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

extern const uint16_t PWM_IMPULSES [];
//= {0, 40, 100, 150, 200, 250, 300, 350, 400, 500, 999};

/**
 * \def KISS_FFT_FFT is KISS FFT computed
 */
#define KISS_FFT_FFT  0
/**
 * \def KISS_FFT_IFFT is KISS IFFT computed
 */
#define KISS_FFT_IFFT 1
/**
 *
 * @param nfft
 * @param inverse_fft
 */
void init_fft (size_t nfft, int inverse_fft);
/**
 *
 * @param input
 * @param output
 * @param length
 */
void real_fft (kiss_fft_cpx * input, kiss_fft_cpx * output, size_t length);
/**
 *
 * @param x
 * @param h
 * @param output
 * @param xlen
 * @param hlen
 */
void fast_conv_f32(floatp_t * x, floatp_t * h, floatp_t * output, size_t xlen,
        size_t hlen);
/**
 *
 */
void clean_fft (void);

/**
 *
 * @param arg
 * @param power
 * @return
 */
floatp_t fast_pow_of(const floatp_t arg, unsigned int power);
/**
 *
 * @param val
 * @return
 */
floatp_t fast_abs(floatp_t val);

#endif /* INC_UTILS_H_ */
