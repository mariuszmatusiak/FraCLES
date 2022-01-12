/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    types.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   types.h header file with type and structs declarations for FraCLES
 * @date    2019/07/07
 ******************************************************************************/

#ifndef INC_TYPES_H_
#define INC_TYPES_H_

#define OPTIMIZE_1
#ifdef OPTIMIZE_1

#ifdef STM32F7

#define    ARM_MATH_CM7
#define    __FPU_PRESENT 1
#define    __DSP_PRESENT 1
#define    __FPU_USED 1
#define ADC_RESOLUTION_BITS 12
#define ADC_RESOLUTION (1 << ADC_RESOLUTION_BITS)
#define DAC_RESOLUTION_BITS 12
#define DAC_RESOLUTION (1 << DAC_RESOLUTION_BITS)

#elif defined STM32H7

#define    ARM_MATH_CM7
#define    __FPU_PRESENT 1
#define    __DSP_PRESENT 1
#ifndef    __FPU_USED
#define    __FPU_USED 1
#endif
#define ADC_RESOLUTION_BITS 16
#define ADC_RESOLUTION (1 << ADC_RESOLUTION_BITS)
#define DAC_RESOLUTION_BITS 12
#define DAC_RESOLUTION (1 << DAC_RESOLUTION_BITS)

#elif defined STM32L1

#define ARM_MATH_CM3
#define __FPU_PRESENT 0

#endif

#include <arm_math.h>
#endif

#ifdef  FIXED_POINT
typedef q31_t floatp_t;
#define Q_M_N 31
#define Q_N   28
#define Q_M   (Q_M_N - Q_N)
const int q_m_n_factor = (1 << Q_N);
#define FPU_CFG "Q31"

#elif defined DOUBLE_PRECISION || defined STM32H7

/**
 *
 */
typedef double floatp_t;
typedef union {
    double f;
    struct
    {
        // IEEE 64-bit float
        uint64_t mantissa : 52;
        unsigned int exponent : 11;
        unsigned int sign : 1;
    } bytes;
} floatp_u;
#define FPU_CFG "F64"
#elif defined STM32F7

typedef float32_t floatp_t;
typedef union {
    float32_t f;
    struct
    {
        // IEEE 32-bit float
        unsigned int mantissa : 23;
        unsigned int exponent : 8;
        unsigned int sign : 1;
    } bytes;
} floatp_u;
#define FPU_CFG "F32"

#else

typedef float32_t floatp_t;
typedef union {
    float32_t f;
    struct
    {
        // IEEE 32-bit float
        unsigned int mantissa : 23;
        unsigned int exponent : 8;
        unsigned int sign : 1;
    } bytes;
} floatp_u;
#define FPU_CFG "F32"

#endif

// max and min ADC input and DAC output values
#ifndef FIXED_POINT

#define MIN_VOLTAGE_VALUE (floatp_t)0.0
#define MAX_VOLTAGE_VALUE (floatp_t)3.3
#define STARTING_POINT    (floatp_t)0.0
#define MAX_DAC_VALUE     (floatp_t)3.3
#define MIN_DAC_VALUE     (floatp_t)0.0

#else

#define MIN_VOLTAGE_VALUE (floatp_t)0
#define MAX_VOLTAGE_VALUE (floatp_t)885837004// 3.3f * 2^Q_N
#define STARTING_POINT    0.0
#define MAX_DAC_VALUE     (floatp_t)885837004
#define MIN_DAC_VALUE     (floatp_t)0.0

#endif

#define MAX_PID_BUFFER_SIZE 500

#ifdef DAC_RESOLUTION

#define MIN_DAC_OUT_VALUE 0
#define MAX_DAC_OUT_VALUE (DAC_RESOLUTION-1)

#else

#define MIN_DAC_OUT_VALUE 0
#define MAX_DAC_OUT_VALUE ((1 << 12) - 1)

#endif

#ifdef ADC_RESOLUTION

#define MIN_ADC_IN_VALUE 0
#define MAX_ADC_IN_VALUE (ADC_RESOLUTION-1)

#else

#define MIN_ADC_IN_VALUE 0
#define MAX_ADC_IN_VALUE ((1 << 12) - 1)

#endif

#ifndef MAX_NUM_DEN_ORDER
#define MAX_NUM_DEN_ORDER 5
#endif

#define true  1u
#define false 0u

/**
 *
 */
enum DISCRETIZATION_METHOD {
    ZOH,    ///< Zero-order hold - assumes the control inputs are piecewise
            ///<                   constant over the sample time Ts.
    FOH,    ///< First-order hold (Triangle approximation, modified) - assumes
            ///<                   the control inputs are piecewise linear over
            ///<                   the sample time T
    IMPULSE,///< Impulse invariant discretization
    TUSTIN, ///< Bilinear (Tustin) method
    MATCHED ///< Zero-pole matching method.
};

/**
 * \brief Struct for 2nd order plant coefficients
 *
 * Realizes the formula:
 *                Kp
 * G(s) = ----------------- * exp(-Td*s)
 *        (1+Tp1*s)(1+Tp2*s)
 *
 * Example parameters from MATLAB:
 *  - Kp  = 1
 *  - Tp1 = 0.31944
 *  - Tp2 = 0.31944
 *  - Td  = 1.0642
 */
typedef struct second_order_dynamic_plant_coefficients {
    float Kp;
    float T_p1;
    float T_p2;
    float T_d;
} c2ndOrderPlantCoeffs;

/**
 *  \brief Struct for coeffs of a discrete plant
 *
 *  Realizes the formula:
 *               b[0]*z^(-nb[0]) + b[1]*z^-nb[1] + b[2]*z^-nb[2]
 *  z^(-delay) * ------------------------------------------------
 *               a[0]*z^(-na[0]) - a[1]*z^-na[1] + a[2]*z^-na[2]
 */
typedef struct f_discrete_dynamic_plant_coefficients {
    //TODO dynamic initialization
    floatp_t  b[MAX_NUM_DEN_ORDER];
    floatp_t  a[MAX_NUM_DEN_ORDER];
    int      nb[MAX_NUM_DEN_ORDER];
    int      na[MAX_NUM_DEN_ORDER];
    int      num_order; //b
    int      den_order; //a
    int      delay;
} discretePlantCoeffs;

/**
 * \brief Structure for operations on fractional-order derivatives and integrals
 *
 * Realizes the formula:
 *  \f$ a_{j} \cdot t_{0} D_{t_{a}}^{a_{n}} \f$
 *
 */
typedef struct fo_differintegral {
    floatp_t a_j; ///< \f$ a_{j} \f$ coefficient
    floatp_t a_n; ///< \f$ a_{n} \f$ order
    floatp_t   h; ///< \f$ h \f$ sample time
    size_t     N; ///< \f$ N \f$ maximum number of samples to compute,
                  ///< \f$ \max{N} = \[\frac{(t_a - t_0)}{h}\] \f$
    // private variables
    // floatp_t m_t_0; // starting time t0 not supported
    // floatp_t m_t_a; // end time      ta not supported
    floatp_t m_h_pow_a_n; ///< Value of \f$ h^{a}_{n} \f$. This parameter is
                          ///< used to optimize computations.
    floatp_t * a; ///< pointer to \f$ a^{\nu}\f$ coefficients array
} fo_differintegral;

typedef struct fotf {
    fo_differintegral b_j[MAX_NUM_DEN_ORDER]; ///< \f$ b_{j} \f$ numerator coeff
                                              ///< of a fractional-order
                                              ///< transfer function
    fo_differintegral a_i[MAX_NUM_DEN_ORDER]; ///< \f$ a_{i} \f$ denominator
                                              ///< coefficients of a fractional-
                                              ///< order transfer function
    size_t num_order_length;    ///< \f$ n \leq N-1 \f$ numerator length
    size_t den_order_length;    ///< \f$ n \leq N-1 \f$ denominator length
    size_t delay;               ///< \f$ e^{-Td} \f$ delay template
} fotf;

typedef struct oustaloupApproximation_Properties
{
     float omega_h; ///< unsigned float upper frequency bandwidth limit
     float omega_l; ///< unsigned float lower frequency bandwidth limit
     float q_order; ///< float fractional order
     uint32_t N;    ///< order of Oustaloup approximation
} oustaloupApproximation_Properties;

/**
 *
 */
typedef struct pid_gains_t {
    floatp_t Kp; ///< \f$ K_P \f$ gain of (VFO)PID
    floatp_t Ki; ///< \f$ K_I \f$ gain of (VFO)PID
    floatp_t Kd; ///< \f$ K_D \f$ gain of (VFO)PID
} pid_gains_t;

typedef struct orders_pair_t {
    floatp_t mu; ///< \f$ \mu \f$ order of (V)FOPID
    floatp_t nu; ///< \f$ \nu \f$ order of (V)FOPID
} orders_pair_t;

typedef struct fopid_parameters_t {
    pid_gains_t   gains; ///< \f$ K_P, K_I, K_D \f$
    orders_pair_t orders;///< \f$ \mu, \nu \f$
} fopid_parameters_t;

typedef enum {
    RESULT_ERROR = -1, ///< Result error code -1
    RESULT_SUCCESS,    ///< Result OK
    RESULT_ERROR_2     ///< Result error code 1
} result_t;

typedef uint8_t bool;

// Typedefs
typedef enum {
    BYPASS_OFF,                 //!< BYPASS, ADC->DAC
    DAC_FIXED_VALUE,            //!< Fixed value is set on MCU DAC output
    CONTROLLER_ENCODER_JUMP,    //!< MCU acts as a controller with jump to fixed
                                //!< encoder position
    CONTROLLER_PID,             //!< MCU acts as a PID controller
    CONTROLLER_FOPID,           //!< MCU acts as a FOPID controller
    CONTROLLER_VFOPID,          //!< MCU acts as a VFOPID controller
    CONTROLLER_VFOPID_2ND_TUNE, //!< MCU acts as a VFOPID controller with an
                                //!< alternative tuning method
    PLANT_INTEGER_ORDER,        //!< MCU acts as integer-order plant
    PLANT_FRACTIONAL_ORDER,     //!< MCU acts as fractional-order plant
} OutputMode;

typedef enum {
    TRAPEZOIDAL,
    SQUARE
} type_integral_t;

extern floatp_t * outputs;
extern floatp_t * inputs;
extern floatp_t * errors;
extern floatp_t * temp_samples;

#endif /* INC_TYPES_H_ */
