/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    dynamic.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   dynamic.h header file with PID and (V)FOPID functions declarations
 * @date    2018/09/24
 ******************************************************************************/

#ifndef DYNAMIC_H_
#define DYNAMIC_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "types.h"
#include "buffers.h"
#include "back_diff_lib.h"

/**
 * The x coefficient of the first_order_response function, set based on the
 *  motor parameters.
 */
#define FIRST_ORDER_X_COEFF 0.1223f

/**
 * The y coefficient of the first_order_response function, set based on the
 * motor parameters.
 */
#define FIRST_ORDER_Y_COEFF 0.8752f

volatile int   xWrIdx;
volatile int   yWrIdx;
volatile int   xRdIdx;
volatile int   yRdIdx;
volatile int   eWrIdx;
size_t    outputs_length;
size_t     inputs_length;
size_t     errors_length;

volatile int m_delayBufferWrIdx;
volatile int m_delayBufferRdIdx;
floatp_t * m_delayBuffer;
size_t  m_delayBufferLength;

fo_differintegral new_fo_differintegral (
        floatp_t a_j, ///< coefficient
        floatp_t a_n, ///< order
        floatp_t   h, ///< sample time
        size_t     N  ///< maximum number of samples to compute,
                      ///< max{N} = [(t_a - t_0) / h]
);

bool new_fotf(
        fotf * empty_fotf,
        size_t num_order_length,
        size_t den_order_length,
        floatp_t * b,
        floatp_t * n_b,
        floatp_t * a,
        floatp_t * n_a,
        size_t delay,
        floatp_t h,
        size_t N_b,
        size_t N_a
);

int  initializeDelayBuffer(size_t length);
void initializeInputOutputBuffers(size_t input_length, size_t output_length);
void initializeErrorBuffer(size_t input_length);
void cleanInputOutputBuffers(void);
void cleanDelayBuffer(void);
void clear_fo_differintegral (fo_differintegral fo);
void clearErrorBuffer(void);
void update_vfopid (fotf * fopid, fopid_parameters_t * new_gains_and_orders,
        floatp_t * a_mu_coeffs, floatp_t * a_nu_coeffs);
int  addNewInputValue(floatp_t value);
int  addNewErrorValue(floatp_t value, bool compensate_error);

floatp_t convertADCToVolts(uint32_t adc_value);
uint32_t convertVoltsToADC(floatp_t volts_value);
floatp_t filter(floatp_t * h, size_t h_length);
floatp_t fotf_lsim(floatp_t y_output[], floatp_t x_input[], fotf * g);

/**
 * max_delay_of_x = nb+delay x[n-N]
 * max_delay_of_y = na y[n-M]
 *        xRd               xWr
 * x[n] = [0  0  0  0  0 ... 1]
 *        yRd    yWr
 * y[n] = [0  0  y1 y2 y3 ... ]
 * @param inputs
 * @param outputs
 * @param plant
 * @return
 */
floatp_t io_plant_response(floatp_t * inputs, floatp_t * outputs,
        const discretePlantCoeffs  * plant);
floatp_t fo_plant_response(floatp_t * inputs, floatp_t * outputs,
        const fotf * plant);
floatp_t fo_plant_response_circ_buff(floatp_t * inputs, floatp_t * outputs,
        const fotf * plant);
floatp_t f_conv(floatp_t * x, floatp_t * h, floatp_t * output, size_t x_size,
        size_t y_size);

/**
 * Calculates first order response based on the previous and current input and
 * output values.
 * @param x_inputs previous and current input values
 * @param y_outputs previous output values
 * @param index current tables index
 * @param totalSize total size of the tables
 * @return new output value
 */
floatp_t first_order_response(floatp_t x_input, floatp_t * y_outputs, const int
        y_index, uint32_t totalSize);
#endif /* DYNAMIC_H_ */
