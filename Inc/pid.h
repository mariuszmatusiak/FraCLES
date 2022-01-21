/*******************************************************************************
 * Copyright (C) 2020 Mariusz Matusiak
 * @file    pid.h - header file
 * @author  Mariusz Matusiak <mmatusiak@iis.p.lodz.pl>
 *                                     ORCID: 0000-0003-2407-3135
 * @brief   pid.h header file with PID and (V)FOPID functions declarations
 * @date    2018/09/24
 ******************************************************************************/

#ifndef PID_H_
#define PID_H_

#include <stdint.h>
#include "types.h"

// Below values taken from Matlab simulations on first_second_order.m

//#define KP_FOR_FIRST_ORDER_ENGINE (1.5f)
#define KP_FOR_FIRST_ORDER_ENGINE (2.0f)
//#define TI_FOR_FIRST_ORDER_ENGINE (0.09f)
#define TI_FOR_FIRST_ORDER_ENGINE (0.028f)
//#define TD_FOR_FIRST_ORDER_ENGINE (0.01f)
#define TD_FOR_FIRST_ORDER_ENGINE (0.007f)

#define ERR_THRESHOLD (0.005)

#define K_ANTI_WINDUP (floatp_t)1

// Define whether apply output values limitations
typedef enum  {
    NONE,
    CONDITIONAL_INTEGRATION_CLAMPING,
    BACK_CALCULATION,
    PADULA_METHOD
} PID_antiwindup_strategy;

typedef struct PID_value {
    floatp_t pid_value;             ///< \f$ u_{min} \leq u(t) \leq u_{max} \f$
    floatp_t p_component;           ///< \f$ u_{P}(t) \f$
    floatp_t i_component;           ///< \f$ u_{I}(t) \f$
    floatp_t d_component;           ///< \f$ u_{D}(t) \f$
    floatp_t pid_not_limited;       ///< \f$ u(t) \f$
    bool isFOPID;                   ///< Determines whether FOPID is still
} PID_value;

typedef struct PID_options {
    bool is_saturation_enabled;     ///< Determines whether controller output
                                    ///< is saturated.
    bool is_i_term_enabled;         ///< Determines whether I term exists.
    bool is_d_term_enabled;         ///< Determines whether D term exists.
    PID_antiwindup_strategy antiwindup_method; ///< Determines whether
                                    ///< back-calculation or other type is used
                                    ///< anti-windup is used.
    floatp_t saturation_lower_limit;    ///< Saturation lower boundary.
    floatp_t saturation_upper_limit;    ///< Saturation upper boundary.
    floatp_t error_threshold;       ///< \f$ E_{ss} \f$ threshold of error in SS
    floatp_t anti_windup_filter;    ///< Value of the anti-windup filter coeff.
                                    ///<    -T_t - tracking time constant e.g.
                                    ///<           \sqrt(T_iT_d)
                                    ///<    -
    floatp_t t_s_sampling_time;     ///< Sampling time \f$ h \f$ [s].
    size_t steady_state_limit;      ///< Maximum transition length up to State.
    type_integral_t integral_type;  ///< Type of the integral to be used.
    bool is_pid_switching_enabled;  ///< Determines whether (V)FOPID shall
                                    ///< change to PID on steady state detection
} PID_options;

typedef struct Controller_state {
    PID_value last_control_output;  //!< Keeps the information about the \f$
                                    //!< u(kh-h) \f$ to keep the integral value
                                    //!< on (V)FOPID -> PID switching
    floatp_t error_integral_approximation;  //!< Keeps the information about
                                    //!< accumulated integral value for (V)FOPID
                                    //!< -> PID switching
    floatp_t last_error_integral_input;   //!< Kp * e(kh-h) / Ti or filtered by
                                    //!< anti-windup block.
    floatp_t previous_error;        //!< Keeps the information about the e(kh-h)
    floatp_t previous_error_abs;    //!< Keeps the information about the |e(kh-
                                    //!< h)| for steady state detection
    uint32_t steady_state_counter;  //!< Keeps the information about the number
                                    //!< of steps in the range +/-E_C to deter-
                                    //!< -mine whether steady state is reached
    bool isPID;                     //!< keeps the information about the status
                                    //!< of the controller: (V)FOPID or PID
    floatp_t pid_reference_value;   //!< r(t)
} Controller_state;

// Global variables
extern volatile floatp_t pid_reference_value;
extern OutputMode output_mode;

// Declarations


/**
 *
 * @param table
 * @param BUFFER_SIZE
 */
void shiftArrayLeft(float * table, const int BUFFER_SIZE);

/**
 *
 * @param table
 * @param index
 * @param element
 * @param BUFFER_SIZE
 */
void addValueToTheTable(float * table, int * index, float element,
        const int BUFFER_SIZE);

/**
 *
 * @param table
 * @param BUFFER_SIZE
 */
void shiftArrayLeftForPIDType(PID_value * table, const int BUFFER_SIZE);

/**
 *
 * @param table
 * @param index
 * @param element
 * @param BUFFER_SIZE
 */
void addValueToTheTableForPIDType(PID_value * table, int * index,
        PID_value element, const int BUFFER_SIZE);

/**
 * PID algorithm with anti-windup filtering, saturation block, scheme block
 * REQUIRED
 * @param error
 * @param pid_params
 * @param opts
 * @return PID_value
 */
PID_value pid(floatp_t curr_error, floatp_t * pid_params, PID_options * opts);

/**
 *
 * @param errors
 * @param fopid
 * @param opts
 * @return
 */
PID_value fractional_pid(floatp_t * errors, fotf * fopid, PID_options * opts);

/**
 *
 */
void fractional_pid_reset(void);

/**
 *
 * @param previous_error
 * @param new_error
 * @param t_s
 * @param type
 * @return
 */
floatp_t calculate_integral_extension(floatp_t previous_error,
        floatp_t new_error, floatp_t t_s, const type_integral_t type);
/**
 * Calculates integral approximation - trapezoidal
 * @param errors
 * @param size
 * @return
 */
floatp_t calculate_integral_approximation(floatp_t * errors,
        fo_differintegral foi, uint32_t index, floatp_t antiWindup);

/**
 *
 * @param result
 * @param error_integral_input
 * @param opts
 * @param controller_gains
 * @param error
 * @return
 */
bool antiwindup_block(PID_value * result, floatp_t * error_integral_input,
        PID_options * opts, floatp_t Ti, floatp_t Td, floatp_t error);
/**
 *
 * @param curr_error
 * @param E
 * @param N
 * @return
 */
bool steadyStateDetector(floatp_t curr_error, const floatp_t E, const size_t N);

#endif /* PID_H_ */
