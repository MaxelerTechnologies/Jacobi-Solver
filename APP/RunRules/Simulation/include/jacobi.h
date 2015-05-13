/**\file */
#ifndef SLIC_DECLARATIONS_jacobi_H
#define SLIC_DECLARATIONS_jacobi_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define jacobi_C (62)
#define jacobi_maxDimLen (8)
#define jacobi_isSimulation (1)
#define jacobi_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_dim Interface Parameter "dim".
 * \param [in] param_equation_num Interface Parameter "equation_num".
 * \param [in] param_max_iter Interface Parameter "max_iter".
 * \param [in] instream_A Stream "A".
 * \param [in] instream_size_A The size of the stream instream_A in bytes.
 * \param [in] instream_B Stream "B".
 * \param [in] instream_size_B The size of the stream instream_B in bytes.
 * \param [in] instream_DiagA Stream "DiagA".
 * \param [in] instream_size_DiagA The size of the stream instream_DiagA in bytes.
 * \param [in] instream_x_trans_init Stream "x_trans_init".
 * \param [in] instream_size_x_trans_init The size of the stream instream_x_trans_init in bytes.
 * \param [out] outstream_error Stream "error".
 * \param [in] outstream_size_error The size of the stream outstream_error in bytes.
 * \param [out] outstream_result Stream "result".
 * \param [in] outstream_size_result The size of the stream outstream_result in bytes.
 */
void jacobi(
	uint64_t param_dim,
	uint64_t param_equation_num,
	uint64_t param_max_iter,
	const void *instream_A,
	size_t instream_size_A,
	const void *instream_B,
	size_t instream_size_B,
	const void *instream_DiagA,
	size_t instream_size_DiagA,
	const void *instream_x_trans_init,
	size_t instream_size_x_trans_init,
	void *outstream_error,
	size_t outstream_size_error,
	void *outstream_result,
	size_t outstream_size_result);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_dim Interface Parameter "dim".
 * \param [in] param_equation_num Interface Parameter "equation_num".
 * \param [in] param_max_iter Interface Parameter "max_iter".
 * \param [in] instream_A Stream "A".
 * \param [in] instream_size_A The size of the stream instream_A in bytes.
 * \param [in] instream_B Stream "B".
 * \param [in] instream_size_B The size of the stream instream_B in bytes.
 * \param [in] instream_DiagA Stream "DiagA".
 * \param [in] instream_size_DiagA The size of the stream instream_DiagA in bytes.
 * \param [in] instream_x_trans_init Stream "x_trans_init".
 * \param [in] instream_size_x_trans_init The size of the stream instream_x_trans_init in bytes.
 * \param [out] outstream_error Stream "error".
 * \param [in] outstream_size_error The size of the stream outstream_error in bytes.
 * \param [out] outstream_result Stream "result".
 * \param [in] outstream_size_result The size of the stream outstream_result in bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *jacobi_nonblock(
	uint64_t param_dim,
	uint64_t param_equation_num,
	uint64_t param_max_iter,
	const void *instream_A,
	size_t instream_size_A,
	const void *instream_B,
	size_t instream_size_B,
	const void *instream_DiagA,
	size_t instream_size_DiagA,
	const void *instream_x_trans_init,
	size_t instream_size_x_trans_init,
	void *outstream_error,
	size_t outstream_size_error,
	void *outstream_result,
	size_t outstream_size_result);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t param_dim; /**<  [in] Interface Parameter "dim". */
	uint64_t param_equation_num; /**<  [in] Interface Parameter "equation_num". */
	uint64_t param_max_iter; /**<  [in] Interface Parameter "max_iter". */
	const void *instream_A; /**<  [in] Stream "A". */
	size_t instream_size_A; /**<  [in] The size of the stream instream_A in bytes. */
	const void *instream_B; /**<  [in] Stream "B". */
	size_t instream_size_B; /**<  [in] The size of the stream instream_B in bytes. */
	const void *instream_DiagA; /**<  [in] Stream "DiagA". */
	size_t instream_size_DiagA; /**<  [in] The size of the stream instream_DiagA in bytes. */
	const void *instream_x_trans_init; /**<  [in] Stream "x_trans_init". */
	size_t instream_size_x_trans_init; /**<  [in] The size of the stream instream_x_trans_init in bytes. */
	void *outstream_error; /**<  [out] Stream "error". */
	size_t outstream_size_error; /**<  [in] The size of the stream outstream_error in bytes. */
	void *outstream_result; /**<  [out] Stream "result". */
	size_t outstream_size_result; /**<  [in] The size of the stream outstream_result in bytes. */
} jacobi_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void jacobi_run(
	max_engine_t *engine,
	jacobi_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *jacobi_run_nonblock(
	max_engine_t *engine,
	jacobi_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void jacobi_run_group(max_group_t *group, jacobi_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *jacobi_run_group_nonblock(max_group_t *group, jacobi_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void jacobi_run_array(max_engarray_t *engarray, jacobi_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *jacobi_run_array_nonblock(max_engarray_t *engarray, jacobi_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* jacobi_convert(max_file_t *maxfile, jacobi_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* jacobi_init(void);

/* Error handling functions */
int jacobi_has_errors(void);
const char* jacobi_get_errors(void);
void jacobi_clear_errors(void);
/* Free statically allocated maxfile data */
void jacobi_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int jacobi_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int jacobi_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_jacobi_H */

