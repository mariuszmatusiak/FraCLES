/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    back_diff_lib.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   back_diff_lib.h header file with the declarations of the functions
 *          designed for calculating the binomial coefficients, backward
 *          differences, fractional derivatives etc.
 * @date    2018/09/24
 ******************************************************************************/

#ifndef BACK_DIFF_LIB_H_
#define BACK_DIFF_LIB_H_


#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "types.h"
#include "buffers.h"
#include "dynamic.h"

/**    USED     **/

/* Maximum measured number of recurrence calls is 7937
 */
/**
 * Calculates array of a coefficients for a given order v and size i.
 * Maximum measured number of recurrence calls is 7937.
 * @param v fractional order
 * @param i size of an array
 * @param a_params_array pointer to initialized array for coefficients
 * @return an element of an array calculated recursively
 */
floatp_t calculateAParamsRecursively(const floatp_t v, unsigned int i,
        floatp_t * a_params_array);
void calculateCParams(const floatp_t v, unsigned int i,
        floatp_t * c_params_array);
/**
 * Calculates backward difference for n samples of inputs using already defined a_params of v order.
 * @param n_samples number of a coefficients
 * @param a_params a coefficients
 * @param reversed_inputs reversed inputs array so the most recent sample is first sample
 * @return backwardDiff fractional order backward difference
 */
floatp_t calculateBackwardDiffStandard(size_t n_samples, floatp_t * a_params,
        floatp_t * inputs, uint32_t j_start);
floatp_t calculateBackwardDiffStandardCircBuff(size_t n_samples,
        floatp_t * a_params, floatp_t * inputs, uint32_t lastElementIdx,
        uint32_t j_start);
/**
 * Calculates backward derivative for n samples of inputs using Horner method and already defined c_params of v order.
 * @param n_samples number of a coefficients
 * @param a_params a coefficients
 * @param reversed_inputs reversed inputs array so the most recent sample is first sample
 * @param v fractional order
 * @param h sampling time
 * @return backwardDiff fractional order backward difference
 */
float calculateBackwardDiffHorner(uint32_t n_samples, float * c_params,
        float * inputs);
/**
 * Calculates backward derivative for n samples of inputs using already defined a_params of v order.
 * @param n_samples number of a coefficients
 * @param a_params a coefficients
 * @param reversed_inputs reversed inputs array so the most recent sample is first sample
 * @param v fractional order
 * @param h sampling time
 * @return backwardDiff fractional order backward difference
 */
floatp_t calculateBackwardDerivativeStandard(size_t n_samples,
        fo_differintegral * fo_differintegral, floatp_t * inputs,
        uint32_t j_start);
floatp_t calculateBackwardDerivativeStandardCircBuff(size_t n_samples,
        fo_differintegral * fo_differintegral, floatp_t * inputs,
        uint32_t lastElementIdx, uint32_t j_start);

#ifdef OPTIMIZE_1

/////////////////////////////////
/// F32 OPTIMIZED
/////////////////////////////////
/**
 * Optimized recursive algorithm using float32_t
 * @param v_order - const fractional positive/negative order
 * @param maximum_index - const index of the last element in the a_parrams_array (SIZE-1)
 * @param a_params_array - pointer to output a_params_array
 * @return last element of the a_params_array
 */
float32_t calculateAParamsFloatRecursivelyOptimized(const float32_t v_order,
        const uint32_t maximum_index, float32_t * a_params_array);

/**
 * Optimized recursive algorithm using q31_t
 * @param v_order - const fractional positive/negative order
 * @param maximum_index - const index of the last element in the a_parrams_array (SIZE-1)
 * @param a_params_array - pointer to output a_params_array
 * @return last element of the a_params_array
 */
q31_t calculateAParamsQ31RecursivelyOptimized(const q31_t v_order,
        const uint32_t maximum_index, q31_t * a_params_array);

/**
 * Optimized algorithm using float32_t and arm_conv_partial_f32
 * @param samples_size - const size of given arrays
 * @param a_params_array - pointer to input a_params_array
 * @param x_inputs - pointer to x_inputs input array
 * @param y_outputs - pointer to y_outputs output array
 * @return last element of y_outputs array
 */
float32_t calculateBackwardDiffStandardOptimized(const uint32_t samples_size,
        float32_t * a_params_array, float32_t * x_inputs, float32_t * y_outputs);

/**
 * Optimized algorithm using q31_t and arm_conv_partial_fast_q31
 * @param samples_size - const size of given arrays
 * @param a_params_array - pointer to input a_params_array
 * @param x_inputs - pointer to x_inputs input array
 * @param y_outputs - pointer to y_outputs output array
 * @return last element of y_outputs array
 */
q31_t calculateBackwardDiffStandardQ31Optimized(const uint32_t samples_size,
        q31_t * a_params_array, q31_t * x_inputs, q31_t * y_outputs);
void arm_Qm_n_to_float(q31_t * x_input_optimized_q31,
        float32_t * x_input_optimized, uint32_t array_size, uint16_t Qm,
        uint16_t Qn);
void arm_float_to_Qm_n(float32_t * x_input_optimized,
        q31_t * x_input_optimized_q31, uint32_t array_size, uint16_t Qm,
        uint16_t Qn);
void calculateAParamsQ31Optimized(const q31_t v_order, q31_t * a_params_array,
        const uint32_t size);
/**
 * Optimized algorithm using float32_t and arm_conv_partial_f32 and arm_scale_f32
 * @param samples_size - const size of given arrays
 * @param a_params_array - pointer to input a_params_array
 * @param x_inputs - pointer to x_inputs input array
 * @param y_outputs - pointer to y_outputs output array
 * @param denominator - const scale factor in a form of 1/powf(base, power)
 * @return last element of y_outputs array
 */
float32_t calculateBackwardDerivativeStandardOptimized(
        const uint32_t samples_size, float32_t * a_params_array,
        float32_t * x_inputs, float32_t * y_outputs,
        const float32_t denominator);

/**
 * Optimized algorithm using q31_t and arm_conv_partial_fast_q31 and arm_scale_q31
 * @param samples_size - const size of given arrays
 * @param a_params_array - pointer to input a_params_array
 * @param x_inputs - pointer to x_inputs input array
 * @param y_outputs - pointer to y_outputs output array
 * @param denominator - const scale factor in a form of 1/powf(base, power)
 * @return last element of y_outputs array
 */
q31_t calculateBackwardDerivativeStandardQ31Optimized(
        const uint32_t samples_size, q31_t * a_params_array, q31_t * x_inputs,
        q31_t * y_outputs, const q31_t denominator);
/////////////////////////////////
#endif

/** IN PROGRESS **/

/**
 * Calculates the value f_k(mu) = 1/k! x d^k/dw^k ((1-w)/(1+w)) | w=0
 * @param order v > 0, u < 0
 * @param k_size
 * @return f_k(mu)
 */
floatp_t powerSeriesExpansionFkCoeff(const floatp_t order, const unsigned int k_size);
void powerSeriesExpansion(floatp_t samplingTime, floatp_t * buffer, int size,
        float v_order);
/**
 * Calculates series of f_k(mu) for set of k values using f_k(mu) = 1/k! x d^k/dw^k ((1-w)/(1+w)) | w=0
 * @param order
 * @param buffer
 * @param k_size
 */
void powerSeriesExpansionFkCoefficients(const floatp_t order, floatp_t * buffer,
        const unsigned int k_size);
void oustaloupSApprox(floatp_t omega_h, floatp_t omega_l,
        const floatp_t q_order);
static inline float oustaloupApproximation_PolesAlphaCoefficient(
        const oustaloupApproximation_Properties properties) {
    float exponent = properties.q_order / (float) properties.N;
    float base = properties.omega_h / properties.omega_l;
    return powf(base, exponent);
}
static inline float oustaloupApproximation_ZeroEtaCoefficient(
        const oustaloupApproximation_Properties properties) {
    float exponent = 1.0f - (properties.q_order / (float) properties.N);
    float base = properties.omega_h / properties.omega_l;
    return powf(base, exponent);
}
/**
 * Calculates series of frequencies of zeros (omega_zi) and poles (omega_pi)
 * for given bandwidth range [properties.omega_l, properties.omega_h], fractional order
 * properties.q_order and order N of Oustaloup approximation.
 * [1] Oustaloup, A. et al, Frequency-Band Complex Noninteger Differentiator: Characterization and Synthesis
 * [2] Merrikh-Bayat, F. et al, Discrete-time fractional-order PID controller: Definition, tuning, digital realization and some applications
 * [3] Oprzedkiewicz, K. et al, An estimation of accuracy of Oustaloup approximation
 * TODO: Discretization
 * @param omega_zi array of zeros frequencies
 * @param omega_pi array of poles frequencies
 * @param properties properties of the approximation
 */
void oustaloupApproximation_ZerosAndPolesFrequencies(float * omega_zi,
        float * omega_pi, const oustaloupApproximation_Properties properties);
float calculateG_knCoefficient(const unsigned int k_size);

/** DEPRECATED **/
unsigned long long calculateFactorial(unsigned long long n);
long long calculateAParamsIntNumerator(int n, unsigned int i);
float calculateAParamsInt(int n, unsigned int i);
float calculateAParamsFloat(float v, unsigned int i);
/**
 * Calculates response of 2nd order dynamic object
 * @param a1 coefficient a1
 * @param a0 coefficient a2
 * @param b1 coefficient b1
 * @param b0 coefficient b0
 * @param inputs array of input values
 * @param outputs array of previous output values
 * @return response of 2nd order dynamic object
 */
float calculateResponseOf2ndOrderObject(float a1, float a0, float b1, float b0,
        float * inputs, float * outputs, uint32_t k);

/**
 * Calculates product of multiplication of two vertexes
 * @param vertex1
 * @param vertex2
 * @param vectors_len
 * @return product of multiplication
 */
float calculateTwoVertexes(float * vertex1, float * vertex2, size_t vectors_len);
float calculateFractionalOrderDiscreteDifferentiatior(uint16_t k_index,
        uint16_t k_index_0, float * a_params, float v, float * inputs,
        uint16_t size, float h_time);
float calculateFractionalOrderDiscreteDifferentiatorForInertia(uint16_t k_index,
        uint16_t k_index_0, uint16_t y_index, float * a_params, float * inputs,
        float * outputs, float v, uint16_t size, float a1, float b1);
void swapArray(float * array, uint16_t len);
uint16_t calculatePIDOutput(uint16_t * inputs, uint16_t * outputs,
        size_t noOfElements, float Kp, float Ki, float Kd, float h,
        float * aParamsNu, float nu, float * aParamsMu, float mu, int e);
void initializeArrayWithUnitJumpI(uint16_t * array, uint16_t len);
void initializeArrayWithLinFunI(uint16_t * array, uint16_t len);
void initializeArrayWithUnitJumpD(float *array, uint16_t len);
void initializeArrayWithLinFunD(float *array, uint16_t len);
void initializeArrayWithUnitJumpWithOffset(uint32_t * array, uint32_t len,
        uint32_t offset);
void initializeArrayWithUnitDumpWithOffset(float * array, uint32_t len,
        uint32_t offset);
#endif /* BACK_DIFF_LIB_H_ */
